breed [bacteria bacterium]               ;; Defines the bacteria breed
breed [antibodies antibody]              ;; Defines the antibody breed
breed [fdcs fdc]                         ;; Defines the FDC breed
breed [naive-b-cells naive-b-cell]                   ;; Defines the naive-b-cell breed
breed [activated-b-cells activated-b-cell]
breed [gc-b-cells gc-b-cell]
breed [sl-plasma-cells sl-plasma-cell]
breed [ll-plasma-cells ll-plasma-cell]
breed [mem-b-cells mem-b-cell]
breed [breg-cells breg-cell]
breed [tfh-cells tfh-cell]
breed [th0-cells th0-cell]
breed [th1-cells th1-cell]
breed [th2-cells th2-cell]

turtles-own [ in-blood bcr isotype csr-bool time-alive tnfa-threshold breg-threshold cd21-level s1pr1-level s1pr2-level cxcr5-level ccr7-level ebi2r-level pro-breg level-of-activation tnf-a-stimulation exposure-number]
activated-b-cells-own [ response-type ]
mem-b-cells-own [ time-in-follicle ]
antibodies-own [antibody-type]
bacteria-own [ epitope-type num-TI-ep num-TD-ep ]
fdcs-own [presented-antigen time-presenting presented-antigen-type responsiveness]
th0-cells-own [r1 r2 rf th1-activation th2-activation tfh-activation]
tfh-cells-own [bcell-binding-status]
th1-cells-own [bcell-binding-status]
th2-cells-own [bcell-binding-status]
patches-own [ patch-type s1p-level cxcl13-level ccl19-level ebi2-level il2 il4 il6 il10 il12 il15 il21 if-g if-a tnf-a tgf-b]  ;; note ccl19 and ccl25 both are used for b-cell localization to b/t border

globals [
  days-passed         ;; Keeps track of the # of days passed since model setup
  average-cd21-expression
  total-num-of-apoptosed-cells
]

;; Most recent calibration of interface sliders
;; cd21-activation-threshold: -7
;; tnfa-apop-threshold: 275
;; breg-diff-threshold: 181
;; phag-il6-burst: 0

; Called when the "setup" button is clicked. Should be the first action by the user.
to setup
  clear-all

  if RandomRuns? = false [random-seed RandomSeed]

  ; Sets up the world structure (lymph node follicle + surrounding paracortex)
  ask patch 0 0 [ask patches in-radius 200  [set patch-type 1 set pcolor gray ]]                            ; Paracortex (outer) zone
  ask patch 0 0 [ask patches in-radius 49  [set patch-type 0 set pcolor black]]                             ; Follicle (inner) zone
  ask patches with [ (pxcor = -50 or pxcor = -49) and abs(pycor) < 5 ] [ set patch-type 2 set pcolor red ]  ; Exit from follicle

  ; Initializes the FDCs
  create-fdcs 100
  ;ask fdcs [ set shape "square" set color brown move-to one-of patches in-radius 30 with [any? fdcs-here = false and any? fdcs-on neighbors = false]]   ; Ensures FDCs don't spawn ontop of eachother
  ask fdcs [ set shape "square" set color brown ]
  ask fdcs [ move-to one-of patches in-radius 30 with [not any? other fdcs in-radius 3]]   ; Ensures FDCs don't spawn ontop of eachother or too close together


  ; Initializes the starting population of Tfh cells
  create-tfh-cells 50
  ask tfh-cells [ move-to one-of patches with [patch-type = 1] set time-alive -1000 set shape "square" set color cyan  set cxcr5-level 11 set ccr7-level 6 set ebi2r-level 5 set bcell-binding-status false]

  create-th1-cells 10
  ask th1-cells [ move-to one-of patches with [patch-type = 1] set time-alive -1000 set shape "square" set color blue  set cxcr5-level 0 set ccr7-level 2 set ebi2r-level 2 set bcell-binding-status false]
  create-th2-cells 50
  ask th2-cells [ move-to one-of patches with [patch-type = 1] set time-alive -1000 set shape "square" set color blue  set cxcr5-level 0 set ccr7-level 2 set ebi2r-level 2 set bcell-binding-status false]
  ;; Initialize global variables and counters
  set days-passed 0

  reset-ticks
end


; Called every tick
to go
  calculate-incoming-tnfa-il6-level
  ;if ticks > 960 [ ;; start simuluating background inflammation on day 20, or 960 ticks, so that the first inoculation isn't affected by background inflammation
    simulate-background-inflammation
  ;]

  ; Calculates the # of days passed from the # of ticks passed
  set days-passed ticks / 48    ;; 1 tick = 30 minutes, so 48 ticks = 1 day

  spawn-b-cell
  ;spawn-th0-cell

  ; Cytokine release from paracortex
  ask patches with [ patch-type = 1 ] [
   set ccl19-level ccl19-level + 2
   set ebi2-level ebi2-level + 2
  ]

  ; Cytokine release from follicle exit
  ask patches with [ patch-type = 2 ] [
   set s1p-level s1p-level + 2
  ]

  ask fdcs [fdc-function]
  ask naive-b-cells [ naive-b-cell-function ]
  ask activated-b-cells [ activated-b-cell-function ]
  ask gc-b-cells [ gc-b-cell-function]
  ask ll-plasma-cells [ll-plasma-cell-function ]
  ask sl-plasma-cells [sl-plasma-cell-function ]
  ask mem-b-cells [mem-b-cell-function]
  ask antibodies [antibodies-function]
  ask breg-cells [ breg-function ]
  ask tfh-cells [ tfh-cell-function ]
  ask th0-cells [th0-cell-function ]
  ask th1-cells [th1-cell-function ]
  ask th2-cells [th2-cell-function ]
  ask bacteria [ bacteria-function ]



  update-chemokine-gradient

  check-overall-cd21-expression

  ; Automatically inoculates at specified ticks if autoinoculate? switch is ON
  if autoinoculate? [
    if ticks = 480 [    ; At 10 days
      auto-inoculate first-exposure-amt    ; Variable from interface input
    ]
    if ticks = 2880 [   ; At 60 days
      auto-inoculate second-exposure-amt   ; Variable from interface input
    ]
    if ticks = 4800 [   ; Stops the run at 100 days
      ;setup
      stop
    ]
  ]

;  if count bacteria > 0 and ticks mod 40 = 0[  ; used to linearly clear the bacteria over time
;    ask one-of bacteria [ die ]
;  ]

  tick
end


; Spawns naive B-cells at a constant rate
to spawn-b-cell
  if ticks mod 10 = 0 [
    create-naive-b-cells 1 [ set shape "circle" set color white set size 1 setxy 49 0
      set time-alive 0
      set bcr random 30
      set isotype "md"              ;; isotype of "md" is IgM/IgD coexpresion. "d" is IgD, "m" is IgM, "a" is IgA, "g" is IgG, "e" is IgE
      set s1pr1-level 0
      set s1pr2-level 0
      set cxcr5-level 16
      set ccr7-level 0
      set ebi2r-level 0
      set cd21-level 0
      set in-blood false

      set tnfa-threshold (bcell-tnfa-apop-threshold + random 50)
      set breg-threshold (bcell-breg-diff-threshold + random 50)
    ]
  ]
end


; Spawns precursor TH0 cells at a constant rate
to spawn-th0-cell
  if ticks mod 20 = 0 [
    create-th0-cells 1 [ set shape "square" set color yellow
      move-to one-of patches with [patch-type = 1]
      set time-alive 0
      set s1pr1-level 0
      set s1pr2-level 0
      set cxcr5-level 0
      set ccr7-level 6
      set ebi2r-level 6
      set in-blood false
    ]
  ]
end


to fdc-function
  ; Secretes the cytokines below
  set cxcl13-level cxcl13-level + 2
  set il4 il4 + 1
  set il6 il6 + 2
  ;set il6 il6 + 1
  set il15 il15 + 2
  set il12 il12 + 2

  ; If FDC is presenting antigen, determines when it will stop presenting it
  if presented-antigen != 0 [     ; presented-antigen = 0 represents no presented antigen
    set time-presenting time-presenting + 1
  ]
  if time-presenting = 300 [      ; Stops presenting the antigen after 300 ticks / 6 days
    set presented-antigen 0
    set color brown
    set presented-antigen-type 0
    set time-presenting 0
  ]

  calculate-fdc-responsiveness
end

to calculate-fdc-responsiveness
  if responsiveness > 0 [
    set responsiveness responsiveness - 0.01
  ]

end

to antibodies-function
  set time-alive time-alive + 1
  if time-alive > 900 [
    die
  ]
end


to naive-b-cell-function
  set cd21-level 300 - (il6 * 50 + il10 * 0)     ; Calculates CD21 exppression based off of il6 and il10 levels (scaled) in current patch

  if patch-type = 2 [         ; Represents naive b cell exiting the follicle
    die
  ]

  ;; Checks if naive b-cell is in contact with an antigen-bound APC or free-floating antigen
  let apc one-of fdcs-here
  let random-encounter random 100
  let antigen one-of bacteria-here
  if (apc != nobody and [presented-antigen] of apc != 0 and [responsiveness] of apc > random-encounter) or antigen != nobody[

    ; Only activates if the B-cell's cd21 level is above the cd21 threshold (bcell-cd21-activation-threshold is defined in the interface slider)
    ifelse cd21-level > bcell-cd21-activation-threshold [
      ;print word (word breed " activating with cd21-level at: ") cd21-level
      set breed activated-b-cells
      set il6 il6 + phag-il6-burst
      set pro-breg 0
      set shape "circle"
      set size 1
      set color yellow
      set csr-bool false
      set time-alive 0
      ifelse antigen != nobody [
        ; Randomly picks whether the B-cell binds to a TI epitope or a TD epitope
        ;let rTI random [num-TI-ep] of antigen
        ;let rTD random [num-TD-ep] of antigen
        let rTI random 100
        let rTD random 100
        ifelse rTI > rTD [
          set response-type 1   ; Response-type of 1 represents a TI response
        ][
          set response-type 2   ; Response-type of 2 represents a TD response

          ; If TD response, B-cell upregulates ccr7 and ebi2r to localize to paracortex
          set ccr7-level 12
          set ebi2r-level 12
        ]
        ask antigen [ die ]
      ][
        if apc != nobody [
          let rTI random 100
          let rTD random 100
          ifelse rTI > rTD [
            set response-type 1
          ][
            set response-type 2
            set ccr7-level 12
            set ebi2r-level 12
          ]
          ;set response-type [presented-antigen-type] of apc
        ]
      ]
    ][
      ;print  word (word breed "NOT activating with cd21-level at: " ) cd21-level
    ]
  ]

  chemotaxis
  move

  check-breg-status   ; Checks level of stimulation of b-reg differentiation
  check-tnf-status    ; Checks level of TNF-a stimulation for apoptosis

  if time-alive > 300 [
    set s1pr1-level s1pr1-level + 0.5 ;; this slowly increases the # of s1p receptors (s1pr) in the naive b cell when the b-cell is old enough
  ]

  set time-alive time-alive + 1
  if time-alive > 1000 [
    die
  ]
end



to breg-function
  ; Secretes the below cytokines
  set il10 il10 + 5
  set tgf-b tgf-b + 1

  chemotaxis
  ;breg-move
  move
  check-tnf-status

  set time-alive time-alive + 1
  if time-alive > 300 [ ; used to be 300
    die
  ]
end


to activated-b-cell-function

  if in-blood = false [
    ; Only performs the following commands if currently inside the follicle, and NOT in the blood/lymph
    if patch-type = 2 [   ; If the cell reaches the follicle exit (patch-type 2)
      set in-blood true
      hide-turtle
    ]

    isotype-switch   ; Determines which isotype to switch to

    ifelse response-type = 2 [
      td-response
    ][
      ifelse response-type = 1 [
        ti-response
      ][
        ; activated-mem-response    ;THIS IS NOT CURRENTLY BEING USED, KEPT IN CASE WE WANT EFFECTOR MEM FUNCTION TO BE DIFFERENT FROM NORMAL B-CELLs
      ]
    ]

    ;check-breg-status

    chemotaxis
    move
  ]

  check-tnf-status

  set time-alive time-alive + 1
  if time-alive > 300 [
    die
  ]
end


; Determines which isotype an activated B-cell will switch to
to isotype-switch
  if csr-bool = false [
      let igM-bucket 0
      let igD-bucket 0
      let igA-bucket 0
      let igG-bucket 0
      let igE-bucket 0

      set igM-bucket il12 + il15 + il6
      ;set igD-bucket   ;seems igD differentiation isnt stimulated by anything
      set igA-bucket il10 + il15 + il21 + tgf-b
      set igG-bucket il4 + il10 + il15 + il21
      set igE-bucket il4 - il12 - if-a - if-g - tgf-b + il21

      let max_index 0
      let mylist (list 3 igM-bucket igA-bucket igG-bucket igE-bucket )
      foreach (list 1 2 3 4) [
        [x] ->
        let value item x mylist
        if value > item max_index mylist [
          set max_index x
        ]
      ]

      if max_index = 1 [
        set csr-bool true
        set isotype "m"
      ]
      if max_index = 2 [
        set csr-bool true
        set isotype "a"
      ]
      if max_index = 3 [
        set csr-bool true
        set isotype "g"
      ]
      if max_index = 4 [
        set csr-bool true
        set isotype "e"
      ]
    ]
end


; T-dependent B-cell response
to td-response
  ; The activated B-cell can get help from either a Tfh cell or a Th2 cell
  let tfh one-of tfh-cells-here
  let th2 one-of th2-cells-here
  ifelse tfh != nobody [
    set breed gc-b-cells
    set pro-breg 0
    set color orange
    set shape "circle"
    set size 1
    set time-alive 0
    create-link-with tfh [ tie ]
    ask tfh [ set ebi2r-level 0 set ccr7-level 0 set bcell-binding-status true]
  ][
   if th2 != nobody [
      set breed gc-b-cells
      set pro-breg 0
      set color orange
      set shape "circle"
      set size 1
      set time-alive 0
      create-link-with th2 [ tie ]
      ask th2 [ set ebi2r-level 0 set ccr7-level 0 set bcell-binding-status true]
    ]
  ]
end

to ti-response
  ; Activated B-cells undergoing TI response secrete TNF-a
  set tnf-a tnf-a + 1

  if time-alive mod 70 = 0 [      ; This is used to represent the rate of proliferation
    ;let proPC (il21 + il10 + if-a + if-g ) * 5    ; Scaled to create reasonable plasma cell populations
    let proPC (il21 + il10 + if-a + if-g ) * 10
    ;let proMem (il21 + il4); * 100
    let proMem 0
    ifelse proPC > proMem [
      hatch-sl-plasma-cells 1 [ set time-alive 0 set color lime + 3 set shape "circle" set size 1 set s1pr1-level 0 set pro-breg 0]
    ][
      ifelse ticks < 2800 [
        hatch-mem-b-cells 1 [set exposure-number 1 set time-alive 0 set color white set shape "target" set s1pr1-level 10 set pro-breg 0 set cd21-level 0 set cxcr5-level 10 ]
        ;hatch-mem-b-cells 1 [set exposure-number 1 set time-alive 0 set color white set shape "target" set s1pr1-level 0 set pro-breg 0 set cd21-level 0 set cxcr5-level 10]
      ][
        hatch-mem-b-cells 1 [set exposure-number 2 set time-alive 0 set color white set shape "target" set s1pr1-level 10 set pro-breg 0 set cd21-level 0 set cxcr5-level 10 ]
        ;hatch-mem-b-cells 1 [set exposure-number 2 set time-alive 0 set color white set shape "target" set s1pr1-level 0 set pro-breg 0 set cd21-level 0 set cxcr5-level 10]
      ]
    ]
  ]
end


; THIS FUNCTION IS NOT CURRENTLY BEING USED
to activated-mem-response
;  set tnf-a tnf-a + 1
;  if counter mod 100 = 0 [
;    hatch-sl-plasma-cells 1 [ set time-alive 0 set color lime + 3 set shape "circle" set size 1 set s1pr1-level 0 set pro-breg 0]
;  ]
end


to gc-b-cell-function
  if in-blood = false [
    if patch-type = 2 [
      set in-blood true
      hide-turtle
    ]

    ; Downregulates ebi2r and ccr7 so it can localize to follicle center again
    set ebi2r-level 0
    set ccr7-level 0

    ; The GC B-cell only moves if it has not yet reached the follicle center
    ifelse distance patch 0 0 > 15 [
      chemotaxis
      gc-move
    ][
      ; Once the GC B-cell reaches the follicle center, it performs the below:
      ;let proPC (il21 + il10 + if-a + if-g ) * 8  ; Scaled for reasonable plasma cell populations
      let proPC2 (il21 + il10 * 2 + if-a + if-g)  ;* 6
      let proMem2 (il21 + il4)
      let proPC random proPC2
      let proMem random proMem2
      ;let proPC random 100        ;randomly choosing if it hatches a llpc or a mem b cell
      ;let proMem random 100

      if time-alive mod 80 = 0 [  ; For now, hard-coded to represent rate of proliferatoin. Should be influenced by cytokines from state diagram
        ifelse proPC > proMem [
          ifelse ticks < 2800 [
            hatch-ll-plasma-cells 1 [ set exposure-number 1 set time-alive 0 set color lime set shape "circle" set size 1 set s1pr1-level 40 set pro-breg 0 set tnfa-threshold tnfa-threshold - 130]
          ] [
            hatch-ll-plasma-cells 1 [ set exposure-number 2 set time-alive 0 set color lime set shape "circle" set size 1 set s1pr1-level 40 set pro-breg 0 set tnfa-threshold tnfa-threshold - 130]
          ]
        ][
          ifelse ticks < 2800 [
            hatch-mem-b-cells 1 [set exposure-number 1 set time-alive 0 set color white set shape "target" set s1pr1-level 10 set pro-breg 0 set cd21-level 0 set cxcr5-level 10 ]
            ;hatch-mem-b-cells 1 [set exposure-number 1 set time-alive 0 set color white set shape "target" set s1pr1-level 0 set pro-breg 0 set cd21-level 0 set cxcr5-level 10]
          ][
            hatch-mem-b-cells 1 [set exposure-number 2 set time-alive 0 set color white set shape "target" set s1pr1-level 10 set pro-breg 0 set cd21-level 0 set cxcr5-level 10 ]
            ;hatch-mem-b-cells 1 [set exposure-number 2 set time-alive 0 set color white set shape "target" set s1pr1-level 0 set pro-breg 0 set cd21-level 0 set cxcr5-level 10]
          ]
        ]
      ]
    ]

    ;check-breg-status
  ]

  check-tnf-status

  set time-alive time-alive + 1
  if time-alive > 700 [
    ask link-neighbors [ set bcell-binding-status false ]  ; Once the GC B-cell dies, frees the bound Tfh or TH2 cell
    die
  ]
end


to sl-plasma-cell-function
  if in-blood = false [
    if patch-type = 2 [
      set in-blood true
      hide-turtle
    ]
    check-breg-status
    chemotaxis
    move
  ]

  if time-alive mod 50 = 0 [  ; For now, hard-coded to represent rate of proliferatoin. Should be influenced by cytokines from state diagram
    hatch-antibodies 1 [ set time-alive 0 set antibody-type isotype set hidden? true ]
  ]

  check-tnf-status

  set time-alive time-alive + 1
  if time-alive > 240 + (il6 + il21) * 10 [
      die
  ]
end


to ll-plasma-cell-function
  ifelse in-blood = false [
    if patch-type = 2 [
      set in-blood true
      hide-turtle
    ]

    check-breg-status

    chemotaxis
    move
  ][
    if time-alive mod 200 = 0 [ ; For now, hard-coded to represent rate of proliferatoin. Should be influenced by cytokines from state diagram
      hatch-antibodies 1 [ set time-alive 0 set antibody-type isotype set hidden? true  ]
    ]
  ]

  check-tnf-status

  set time-alive time-alive + 1
  if time-alive > 8000 + (il6 + il21) * 10 [
      die
  ]


end


; Memory B-cell function is set up to be very similar to naive-B-cell function
to mem-b-cell-function
  ; Sets the cd21 expression level based off of il6 and il10 in current patch
  set cd21-level 300 - (il6 * 50 + il10 * 0)

  if patch-type = 2 [
    setxy 49 0
    set pro-breg 0
    set s1pr1-level 0
    set time-in-follicle 0
  ]

  set time-in-follicle time-in-follicle + 1
  if time-in-follicle > 300 [
    set s1pr1-level s1pr1-level + 0.5 ;; this slowly increases the # of s1p receptors (s1pr) in the naive b cell when the b-cell is old enough
  ]

  set time-alive time-alive + 1
  if time-alive > 15000 [
    die
  ]

  let apc one-of fdcs-here
  let random-encounter random 100
  let antigen one-of bacteria-here
  if (apc != nobody and [presented-antigen] of apc != 0 and [responsiveness] of apc > random-encounter) or antigen != nobody[
    ifelse cd21-level > bcell-cd21-activation-threshold [
      ;print word (word breed " activating with cd21-level at: ") cd21-level
      set breed activated-b-cells
      set il6 il6 + phag-il6-burst
      set cxcr5-level 16
      set s1pr1-level 0
      set pro-breg 0
      set shape "circle"
      set size 1
      set color yellow
      set csr-bool false
      set time-alive 100
      ifelse antigen != nobody [
        ;let rTI random [num-TI-ep] of antigen
        ;let rTD random [num-TD-ep] of antigen
        let rTI random 100
        let rTD random 100
        ifelse rTI > rTD [
          set response-type 1   ; 1 is TI response
        ][
          set response-type 2   ; 2 is TD
          set ccr7-level 12
          set ebi2r-level 12
        ]
        ask antigen [ die ]
      ][
        if apc != nobody [
          let rTI random 100
          let rTD random 100
          ifelse rTI > rTD [
            set response-type 1
          ][
            set response-type 2
            set ccr7-level 12
            set ebi2r-level 12
          ]
          ;set response-type [presented-antigen-type] of apc
        ]
      ]
    ][
      ;print  word  (word breed " NOT activating with cd21-level at: ") cd21-level
    ]
  ]
  check-breg-status

  chemotaxis
  move

  check-tnf-status


end


to th0-cell-function
  let pro-TH1 (il12 + if-g) * 100
  let pro-TH2 (il10 + il4) * 100
  let pro-TFH (il21 + il12) * 100
  let rTH1 random pro-TH1
  let rTH2 random pro-TH2
  let rTFH random pro-TFH
  set r1 rTH1
  set r2 rTH2
  set rf rTFH
  if rTH1 > rTH2 and rTH1 > rTFH [
    set th1-activation th1-activation + 1
  ]
  if rTH2 > rTH1 and rTH2 > rTFH [
    set th2-activation th2-activation + 1
  ]
  if rTFH > rTH1 and rTFH > rTH2 [
    set tfh-activation tfh-activation + 1
  ]

  ifelse th1-activation >= 20 [
    set breed TH1-cells
    set color blue
    set time-alive 0
    set size 1
    set shape "circle"
  ][
    ifelse th2-activation >= 20 [
      set breed th2-cells
      set color blue
      set size 1
      set shape "circle"
      set time-alive 0
      set bcell-binding-status false
    ][
      if tfh-activation >= 20 [
        set breed tfh-cells
        set cxcr5-level 10
        set color cyan
        set shape "circle"
        set size 1
        set time-alive 0
        set bcell-binding-status false
        set cxcr5-level 10
      ]
    ]
  ]

  chemotaxis
  move

  set time-alive time-alive + 1
  if time-alive > 300
    [die]

end

to tfh-cell-function
  if distance patch 0 0 > 20 or bcell-binding-status = false [
    chemotaxis
    move
  ]

  ; Secretes the following cytokines
  set il21 il21 + 1
  set il4 il4 + 1
  set il2 il2 + 1
  set il10 il10 + 1

  set time-alive time-alive + 1
  ;if time-alive > 500
    ;[die]
end


to th1-cell-function
  chemotaxis
  move

  ; Secretes the following cytokines
  set if-g if-g + 1

  set time-alive time-alive + 1
  ;if time-alive > 500
    ;[die]
end

to th2-cell-function
  if distance patch 0 0 > 20 or bcell-binding-status = false [
    chemotaxis
    move
  ]

  ; Secretes the following cytokines
  set il4 il4 + 1
  set il10 il10 + 1

  set time-alive time-alive + 1
  ;if time-alive > 500 [
    ;die
  ;]
end


to bacteria-function
  if patch-type = 2 [ ;; for bacteria, im having them recirculate through blood. when recirculating, they can either just go back into LN, or can be captured by FDC. random chance of either
    ; Randomly chooses if bacteria will recirculate or will be captured by FDC
    let x random 2
    ifelse x = 0 [
      setxy 49 0
      set s1pr1-level 0
    ][
      if any? fdcs with [presented-antigen = 0] [
        ask one-of fdcs with [presented-antigen = 0] [
          set time-presenting 0
          set presented-antigen bacteria-epitope-type
          set color red
          let rTI random number-of-TI-epitopes
          let rTD random number-of-TD-epitopes
          ifelse rTI > rTD [
            set presented-antigen-type 1   ;; 1 is TI epitope to be presented
          ][
            set presented-antigen-type 2    ;; 2 is TD epitope to be presented
          ]
        ]
        die
      ]
    ]
  ]
  set s1pr1-level s1pr1-level + 0.004
  chemotaxis
  move
end


; Checks level of stimulation to differentiate into a regulatory B-cell
to check-breg-status
  ifelse pro-breg > breg-threshold [
    ;print word breed " turned into breg"
    set breed breg-cells
    set size 1
    set shape "circle"
    set color violet
    set s1pr1-level 0 ; to keep breg cell localized to where all the action is
    set cxcr5-level 10
    set time-alive 0

  ][
    set pro-breg (il6 + il21) * 45
  ]
end


; Checks level of tnf-a stimulation for apoptosis
to check-tnf-status
  set tnf-a tnf-a - 0.01 ;;represents consumption of TNFa
  set tnf-a-stimulation 100 * tnf-a

  if tnf-a-stimulation > tnfa-threshold [
    ;print word breed "APOPTOSE with tnf: " print tnf-a-stimulation
    set total-num-of-apoptosed-cells total-num-of-apoptosed-cells + 1
    die
  ]
end


; Diffuses and evaporates all the chemokines
to update-chemokine-gradient
  diffuse cxcl13-level 1   ;; determines the mobility/solubility of cxcl13
  diffuse ccl19-level 1
  diffuse s1p-level 1
  diffuse ebi2-level 1
  diffuse il2 1
  diffuse il4 1
  diffuse il6 1
  diffuse il10 1
  diffuse il12 1
  diffuse il15 1
  diffuse il21 1
  diffuse if-g 1
  diffuse if-a 1
  diffuse tnf-a 1
  diffuse tgf-b 1

  ask patches [
    set cxcl13-level cxcl13-level * 0.9  ;; takes into account protease-driven degradation of cxcl13
    set ccl19-level ccl19-level * 0.9
    set ebi2-level ebi2-level * 0.9
    set s1p-level s1p-level * 0.9
    set il2 il2 * 0.9
    set il4 il4 * 0.9
    set il6 il6 * 0.95
    set il10 il10 * 0.95
    set il12 il12 * 0.9
    set il15 il15 * 0.9
    set il21 il21 * 0.9
    set if-g if-g * 0.9
    set if-a if-a * 0.9
    set tnf-a tnf-a * 0.9
    set tgf-b tgf-b * 0.9

    ; Used to visualize the cytokine levels in the follicle, based off of the interface chooser
    if patch-type = 0 [
      ;let total-cytokine-level il2 + il4 + il6 + il10 + il12 + il15 + il21 + tnf-a + tgf-b + if-a + if-g
      if cytokine-to-visualize = "tnf-a" [
        set pcolor scale-color green tnf-a 0.1 3  ;;used to visualize cxcl13 or ccl19 gradient
      ]
      if cytokine-to-visualize = "il6" [
        set pcolor scale-color green il6 0.1 3  ;;used to visualize cxcl13 or ccl19 gradient
      ]
      if cytokine-to-visualize = "il10" [
        set pcolor scale-color green il10 0.1 3  ;;used to visualize cxcl13 or ccl19 gradient
      ]
      if cytokine-to-visualize = "s1p" [
        set pcolor scale-color green s1p-level 0.01 3  ;;used to visualize cxcl13 or ccl19 gradient
      ]
    ]
  ]
end


;This inoculate function is called when the user clicks the "inoculate" button in the interface. It manually adds bacteria into the system
to inoculate
  ask up-to-n-of (number-of-bacteria / 2) fdcs [
    set time-presenting 0
    set presented-antigen bacteria-epitope-type
   ;set color 15 + (presented-antigen - 1) * 30
    set color red

    let rTI random number-of-TI-epitopes
    let rTD random number-of-TD-epitopes
    ifelse rTI > rTD [
      set presented-antigen-type 1   ;; 1 is TI
    ][
      set presented-antigen-type 2    ;; 2 is TD
    ]
  ]

  create-bacteria (number-of-bacteria / 2) [                            ;; Creates bacteria. "number-of-bacteria" is a variable controlled by an interface slider
    ;set color 15 + (bacteria-epitope-type - 1) * 30               ;; Sets the color of the bacteria based on epitope type. Uses netlogo's 0-139 color scale (integer values)
    set color red
    set shape "bug"
    set size 2
    setxy 49 0
    set s1pr1-level 8
    set time-alive 0
    set in-blood false
    set epitope-type bacteria-epitope-type                        ;; Sets the bacteria's epitope-type. "bacteria-epitope-type" is a value is from an interface slider
    set num-TI-ep number-of-TI-epitopes
    set num-TD-ep number-of-TD-epitopes
  ]
end


; This function is called when the AutoInoculate? interface chooser is ON. It is called in the go function
to auto-inoculate [num-bac]

  ;ask up-to-n-of max (list 60 num-bac) fdcs [
  let num-phagocytosed-bac min (list round (num-bac * 0.8) 100)
  let num-free-bac num-bac - num-phagocytosed-bac

  ask up-to-n-of num-phagocytosed-bac fdcs [

    set responsiveness min (list (responsiveness + 50) 100)
    ;set responsiveness 100
    set time-presenting 0
    set presented-antigen bacteria-epitope-type
   ;set color 15 + (presented-antigen - 1) * 30
    set color red

    let rTI random number-of-TI-epitopes
    let rTD random number-of-TD-epitopes
    ifelse rTI > rTD [
      set presented-antigen-type 1   ;; 1 is TI
    ][
      set presented-antigen-type 2    ;; 2 is TD
    ]
  ]

  create-bacteria num-free-bac [                            ;; Creates bacteria. "number-of-bacteria" is a variable controlled by an interface slider
    ;set color 15 + (bacteria-epitope-type - 1) * 30               ;; Sets the color of the bacteria based on epitope type. Uses netlogo's 0-139 color scale (integer values)
    set color red
    set shape "bug"
    set size 2
    setxy 49 0
    ;set s1pr1-level 8
    set time-alive 0
    set in-blood false
    set epitope-type bacteria-epitope-type                        ;; Sets the bacteria's epitope-type. "bacteria-epitope-type" is a value is from an interface slider
    set num-TI-ep number-of-TI-epitopes
    set num-TD-ep number-of-TD-epitopes
  ]
end


; Calculates the amount of tnf-a and il6 that is evenly distributed in the follicle based off of # of bacteria still present in the system
to calculate-incoming-tnfa-il6-level
  ask patches [set tnf-a tnf-a + ((count bacteria) / 500)]
  ask patches [set il6 il6 + ((count bacteria) / 500)]
end

to simulate-background-inflammation
  ask patches [set tnf-a tnf-a + background-tnfa]
  ask patches [set il6 il6 + background-il6]
end


to check-overall-cd21-expression
  if ticks mod 20 = 0 [     ; Only calculating avg CD21 expression every 50 ticks to increase run speed
    set average-cd21-expression 0
    ask naive-b-cells [
      set average-cd21-expression average-cd21-expression + cd21-level
    ]
    ask mem-b-cells [
      set average-cd21-expression average-cd21-expression + cd21-level
    ]

    ifelse count naive-b-cells != 0 or count mem-b-cells != 0 [
      set average-cd21-expression average-cd21-expression / (count naive-b-cells + count mem-b-cells)
    ][
      set average-cd21-expression 0
    ]
  ]

end


;; Rotates turtle to face direction based off of chemokine gradients
to chemotaxis
  let rt-turn 0
  let lt-turn 0


  let s1pr1-weight s1pr1-level / 100
  let max-s1p-patch max-one-of neighbors [s1p-level]  ;; or neighbors4
  let angle-to-s1p (towards max-s1p-patch)
  let cur-angle heading
  let x angle-to-s1p - cur-angle
  if x < 0 [
    set x x + 360
  ]
  let y 360 - x
  ifelse x < y
  [ set rt-turn rt-turn + x * s1pr1-weight ]
  [ set lt-turn lt-turn + y * s1pr1-weight ]

  let s1pr2-weight s1pr2-level / 100
  let max-s1pr2-patch max-one-of neighbors [s1p-level]  ;; or neighbors4
  let angle-to-s1pr2 (towards max-s1pr2-patch)
  set cur-angle heading
  set x angle-to-s1pr2 - cur-angle
  if x < 0 [
    set x x + 360
  ]
  set y 360 - x
  ifelse x < y
  [ set rt-turn rt-turn + x * s1pr2-weight ]
  [ set lt-turn lt-turn + y * s1pr2-weight ]

  let cxcr5-weight cxcr5-level / 100
  let max-cxcl13-patch max-one-of neighbors [cxcl13-level]  ;; or neighbors4
  let angle-to-cxcl13 (towards max-cxcl13-patch)
  set cur-angle heading
  set x angle-to-cxcl13 - cur-angle
  if x < 0 [
    set x x + 360
  ]
  set y 360 - x
  ifelse x < y
  [ set rt-turn rt-turn + x * cxcr5-weight ]
  [ set lt-turn lt-turn + y * cxcr5-weight ]

  let ccr7-weight ccr7-level / 100
  let max-ccr7-patch max-one-of neighbors [ccl19-level]  ;; or neighbors4
  let angle-to-ccr7 (towards max-ccr7-patch)
  set cur-angle heading
  set x angle-to-ccr7 - cur-angle
  if x < 0 [
    set x x + 360
  ]
  set y 360 - x
  ifelse x < y
  [ set rt-turn rt-turn + x * ccr7-weight ]
  [ set lt-turn lt-turn + y * ccr7-weight ]

  let ebi2r-weight ebi2r-level / 100
  let max-ebi2r-patch max-one-of neighbors [ebi2-level]  ;; or neighbors4
  let angle-to-ebi2r (towards max-ebi2r-patch)
  set cur-angle heading
  set x angle-to-ebi2r - cur-angle
  if x < 0 [
    set x x + 360
  ]
  set y 360 - x
  ifelse x < y
  [ set rt-turn rt-turn + x * ebi2r-weight ]
  [ set lt-turn lt-turn + y * ebi2r-weight ]

  rt rt-turn
  lt lt-turn
end

;; Moves turtle forward one step with a random turn included
to move
  rt random 50
  lt random 50
  fd 1
end


; Similar to move, but just slower. Used to further delay onset/seeding of the GC response
to gc-move
  rt random 50
  lt random 50
  fd 0.5
end

to breg-move
  rt random 50
  lt random 50
  if patch-ahead 1 != nobody [
    if [patch-type] of patch-ahead 1 = 0 [
      fd 1
    ]
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
797
10
1247
461
-1
-1
4.38
1
10
1
1
1
0
0
0
1
-50
50
-50
50
1
1
1
ticks
30.0

BUTTON
13
11
70
49
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
78
11
133
50
NIL
go\n
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
696
29
777
65
Add Antigen
inoculate
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
325
10
503
43
number-of-bacteria
number-of-bacteria
0
300
50.0
1
1
NIL
HORIZONTAL

SLIDER
327
50
504
83
bacteria-epitope-type
bacteria-epitope-type
1
30
30.0
1
1
NIL
HORIZONTAL

SLIDER
514
13
681
46
number-of-TD-epitopes
number-of-TD-epitopes
0
10
10.0
1
1
NIL
HORIZONTAL

SLIDER
513
48
679
81
number-of-TI-epitopes
number-of-TI-epitopes
0
10
10.0
1
1
NIL
HORIZONTAL

PLOT
289
354
518
500
Ag-Specific B-Cell Populations
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"LLPCs" 1.0 0 -14439633 true "" "plot count ll-plasma-cells with [exposure-number = 1]"
"SLPCs" 1.0 0 -8330359 true "" "plot count sl-plasma-cells"
"Mem B-Cells" 1.0 0 -7500403 true "" "plot count mem-b-cells"
"pen-4" 1.0 0 -8630108 true "" "plot count breg-cells"

PLOT
289
86
449
206
antibody
NIL
Ab Level
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count antibodies"

PLOT
624
84
784
204
responsiveness
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ";plot sum ([il4] of patches)"
"pen-1" 1.0 0 -7500403 true "" "plot mean [responsiveness] of fdcs"

CHOOSER
807
481
963
526
cytokine-to-visualize
cytokine-to-visualize
"none" "tnf-a" "il6" "il10" "s1p"
0

MONITOR
145
10
289
55
Number of Days Elapsed
days-passed
1
1
11

SWITCH
17
79
140
112
RandomRuns?
RandomRuns?
0
1
-1000

SLIDER
147
79
264
112
RandomSeed
RandomSeed
0
1000
578.0
1
1
NIL
HORIZONTAL

SLIDER
16
327
277
360
bcell-cd21-activation-threshold
bcell-cd21-activation-threshold
0
200
97.0
1
1
NIL
HORIZONTAL

MONITOR
980
478
1082
523
NIL
count bacteria
17
1
11

SWITCH
17
121
265
154
AutoInoculate?
AutoInoculate?
0
1
-1000

SLIDER
16
369
269
402
bcell-tnfa-apop-threshold
bcell-tnfa-apop-threshold
0
500
275.0
1
1
NIL
HORIZONTAL

SLIDER
17
411
270
444
bcell-breg-diff-threshold
bcell-breg-diff-threshold
0
500
181.0
1
1
NIL
HORIZONTAL

PLOT
452
84
617
206
Total IL-10 Levels
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot sum ([il10] of patches)"

INPUTBOX
17
153
166
213
first-exposure-amt
50.0
1
0
Number

INPUTBOX
138
153
287
213
second-exposure-amt
50.0
1
0
Number

PLOT
14
210
174
330
Average CD21 Expression
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot average-cd21-expression"

PLOT
175
210
338
330
Total Apoptosis
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot total-num-of-apoptosed-cells"

PLOT
287
523
785
673
Total Lymphocyte Counts
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Lymphocytes" 1.0 0 -2674135 true "" "plot (count mem-b-cells) + (count ll-plasma-cells) + (count sl-plasma-cells ) + (count breg-cells) + (count activated-b-cells) + (count gc-b-cells)"

SLIDER
18
457
267
490
phag-il6-burst
phag-il6-burst
0
20
0.0
1
1
NIL
HORIZONTAL

PLOT
595
207
795
343
LLPCs
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count (ll-plasma-cells with [exposure-number = 1])"
"pen-1" 1.0 0 -2674135 true "" "plot count (ll-plasma-cells with [exposure-number = 2])"

PLOT
393
204
593
354
Mem B Cell Response
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count (mem-b-cells with [exposure-number = 1])"
"pen-1" 1.0 0 -2674135 true "" "plot count (mem-b-cells with [exposure-number = 2])"

PLOT
591
342
788
477
gc b cell count
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"pen-0" 1.0 0 -7500403 true "" "plot count gc-b-cells"

PLOT
813
544
1013
694
IL6
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" ";plot (count tfh-cells with [bcell-binding-status = false]) + (count th2-cells with [bcell-binding-status = false])  "
"TI antigen" 1.0 0 -7500403 true "" ";plot count fdcs with [presented-antigen-type = 1]"
"TD antigen" 1.0 0 -2674135 true "" ";plot count fdcs with [presented-antigen-type = 2]"
"pen-3" 1.0 0 -955883 true "" "plot sum [il6] of patches"

INPUTBOX
20
518
126
578
background-tnfa
0.0
1
0
Number

INPUTBOX
141
519
290
579
background-il6
0.0
1
0
Number

@#$#@#$#@
## Description of the model


- The system takes place in a two-dimensional "square" of blood. Therefore, each patch represents an xy location in the blood
- In the blood, there are 3 types of objects: B-cells, bacteria, and antibodies
	- B-cells are circles
	- Bacteria are stars
	- Antibodies are arrowheads
- Each B-cell expresses a single BCR
- Each bacteria is simplified to express only a single epitope
- Each antibody can bind to a single bacterial epitope

 

## How the model works

- All B-cells, bacteria, and antibodies move in random directions in the system to simulate how they might "float" in the blood
- Bacteria behavior:
	- Each bacteria moves in a random location
- B-cell behavior:
	- B-cells move in random directions
	- If it touches a bacteria, it allows its BCR to interact with the bacteria
	- If the bacteria's epitope happens to fit into the BCR, then the bacteria is phagocytosed and killed
	- Then, the B-cell produces antibodies with the same specificity as the bound BCR
- Antibody behavior:
	- Antibodies move in a random direction
	- Similar to a B-cell, if an antibody touches a bacteria, it neutralizes and kills that bacteria


## How to use it

1. Use the starting-num-b-cells slider to set the initial number of B-cells to add to the system
2. Click the setup button to add the B-cells
3. Click the go button
4. Use the bacteria-epitope-type slider to set the single epitope-type expressed by the bacteria you want to inoculate into the system
5. Click the inoculate button to add the bacteria

Notes
- The graph labeled "Number of Live Bacteria" monitors the number of live bacteria in the system at any given point in time.
- The output box labeled "Duration of Infection" measures how long it takes the immune system to clear all the bacteria for any given "inoculation".

## Things to notice

When a bacteria with a given epitope-type is first introduced to the system, it takes a long time for the immune system to clear it. However, as antibodies are synthesized, the second inoculation of that same epitope-type bacteria will result in a shorter time to clear the infection. On the same thread, introducing a different epitope-type bacteria will result in a longer clear time, as no antibodies have been synthesized yet.

## Things to try

- Try inoculating the system with bacteria of any epitope-type value, and then once all the bacteria die, inoculate the system again with the same epitope-type. The time it takes for the infection to clear should be shorter, as shown in the output box labeled "Duration of Infection".
- Then, try inoculating the system with bacteria of a different epitope-type. The duration to clear the infection should be longer.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Parameter Sensitivity Analysis" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="4800"/>
    <metric>average-cd21-expression</metric>
    <metric>total-num-of-apoptosed-cells</metric>
    <metric>sum ([il10] of patches)</metric>
    <metric>count sl-plasma-cells</metric>
    <metric>count (ll-plasma-cells with [exposure-number = 1])</metric>
    <metric>count (ll-plasma-cells with [exposure-number = 2])</metric>
    <metric>count (mem-b-cells with [exposure-number = 1])</metric>
    <metric>count (mem-b-cells with [exposure-number = 2])</metric>
    <metric>count breg-cells</metric>
    <metric>sum ([tnf-a] of patches)</metric>
    <metric>sum ([il6] of patches)</metric>
    <metric>count bacteria</metric>
    <metric>count mem-b-cells + count sl-plasma-cells + count ll-plasma-cells + count breg-cells + count activated-b-cells + count gc-b-cells</metric>
    <enumeratedValueSet variable="number-of-TD-epitopes">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cytokine">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-TI-epitopes">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-bacteria">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bacteria-epitope-type">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RandomRuns?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="bcell-cd21-activation-threshold" first="50" step="10" last="150"/>
    <steppedValueSet variable="bcell-tnfa-apop-threshold" first="240" step="10" last="310"/>
    <steppedValueSet variable="bcell-breg-diff-threshold" first="140" step="10" last="210"/>
    <enumeratedValueSet variable="phag-il6-burst">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="first-exposure-amt">
      <value value="50"/>
    </enumeratedValueSet>
    <steppedValueSet variable="second-exposure-amt" first="50" step="150" last="200"/>
    <enumeratedValueSet variable="autoinoculate?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-tnfa">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-il6">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Calibrated 10 Runs" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="4800"/>
    <metric>average-cd21-expression</metric>
    <metric>total-num-of-apoptosed-cells</metric>
    <metric>sum ([il10] of patches)</metric>
    <metric>count sl-plasma-cells</metric>
    <metric>count (ll-plasma-cells with [exposure-number = 1])</metric>
    <metric>count (ll-plasma-cells with [exposure-number = 2])</metric>
    <metric>count (mem-b-cells with [exposure-number = 1])</metric>
    <metric>count (mem-b-cells with [exposure-number = 2])</metric>
    <metric>count breg-cells</metric>
    <metric>sum ([tnf-a] of patches)</metric>
    <metric>sum ([il6] of patches)</metric>
    <metric>count bacteria</metric>
    <metric>count mem-b-cells + count sl-plasma-cells + count ll-plasma-cells + count breg-cells + count activated-b-cells + count gc-b-cells</metric>
    <enumeratedValueSet variable="first-exposure-amt">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cytokine-to-visualize">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-TI-epitopes">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bcell-cd21-activation-threshold">
      <value value="-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-bacteria">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bcell-breg-diff-threshold">
      <value value="181"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bacteria-epitope-type">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RandomSeed">
      <value value="578"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-TD-epitopes">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phag-il6-burst">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RandomRuns?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AutoInoculate?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bcell-tnfa-apop-threshold">
      <value value="275"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="second-exposure-amt">
      <value value="50"/>
      <value value="200"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Calibrated 10 Runs Sepsis IL6 TNFa Analysis" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="4800"/>
    <metric>sum ([tnf-a] of patches)</metric>
    <metric>sum ([il6] of patches)</metric>
    <metric>count bacteria</metric>
    <enumeratedValueSet variable="first-exposure-amt">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cytokine-to-visualize">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-TI-epitopes">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bcell-cd21-activation-threshold">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-bacteria">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bcell-breg-diff-threshold">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bacteria-epitope-type">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RandomSeed">
      <value value="578"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-TD-epitopes">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phag-il6-burst">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RandomRuns?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AutoInoculate?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bcell-tnfa-apop-threshold">
      <value value="263"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="second-exposure-amt">
      <value value="300"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="LLPC and Mem B Cells 1st/2nd Exposures Separated" repetitions="10" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="4800"/>
    <metric>count (ll-plasma-cells with [exposure-number = 2])</metric>
    <metric>count (ll-plasma-cells with [exposure-number = 1])</metric>
    <metric>count (mem-b-cells with [exposure-number = 2])</metric>
    <metric>count (mem-b-cells with [exposure-number = 1])</metric>
    <enumeratedValueSet variable="first-exposure-amt">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cytokine-to-visualize">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-TI-epitopes">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bcell-cd21-activation-threshold">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-bacteria">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bcell-breg-diff-threshold">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bacteria-epitope-type">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RandomSeed">
      <value value="578"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-TD-epitopes">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phag-il6-burst">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RandomRuns?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AutoInoculate?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bcell-tnfa-apop-threshold">
      <value value="263"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="second-exposure-amt">
      <value value="30"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Background Inflammation - 1st Normal and 2nd Normal Antigen Stimuli" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>average-cd21-expression</metric>
    <metric>total-num-of-apoptosed-cells</metric>
    <metric>sum ([il10] of patches)</metric>
    <metric>count sl-plasma-cells</metric>
    <metric>count (ll-plasma-cells with [exposure-number = 1])</metric>
    <metric>count (ll-plasma-cells with [exposure-number = 2])</metric>
    <metric>count (mem-b-cells with [exposure-number = 1])</metric>
    <metric>count (mem-b-cells with [exposure-number = 2])</metric>
    <metric>count breg-cells</metric>
    <metric>sum ([tnf-a] of patches)</metric>
    <metric>sum ([il6] of patches)</metric>
    <metric>count bacteria</metric>
    <metric>count mem-b-cells + count sl-plasma-cells + count ll-plasma-cells + count breg-cells + count activated-b-cells + count gc-b-cells</metric>
    <enumeratedValueSet variable="first-exposure-amt">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cytokine-to-visualize">
      <value value="&quot;none&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="background-tnfa" first="0" step="0.01" last="0.1"/>
    <enumeratedValueSet variable="number-of-TI-epitopes">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bcell-cd21-activation-threshold">
      <value value="-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="number-of-bacteria">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bcell-breg-diff-threshold">
      <value value="181"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bacteria-epitope-type">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RandomSeed">
      <value value="578"/>
    </enumeratedValueSet>
    <steppedValueSet variable="background-il6" first="0" step="0.01" last="0.1"/>
    <enumeratedValueSet variable="number-of-TD-epitopes">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phag-il6-burst">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AutoInoculate?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RandomRuns?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bcell-tnfa-apop-threshold">
      <value value="275"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="second-exposure-amt">
      <value value="50"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
