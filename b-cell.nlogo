extensions [ vid ]
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


turtles-own [ in-blood bcr isotype csr-bool time-alive cd21-level s1pr1-level s1pr2-level cxcr5-level ccr7-level ebi2r-level pro-breg level-of-activation tnf-a-stimulation]
activated-b-cells-own [ response-type ]
;naive-b-cells-own [cd21-level]
mem-b-cells-own [time-in-blood]
antibodies-own [antibody-type]
bacteria-own [ epitope-type num-TI-ag num-TD-ag ]
fdcs-own [presented-antigen time-presenting presented-antigen-type ]
th0-cells-own [r1 r2 rf th1-activation th2-activation tfh-activation]
tfh-cells-own [bcell-binding-status]
th1-cells-own [bcell-binding-status]
th2-cells-own [bcell-binding-status]
patches-own [ patch-type s1p-level cxcl13-level ccl19-level ebi2-level il2 il4 il6 il10 il12 il15 il21 if-g if-a tnf-a tgf-b]  ;; note ccl19 and ccl25 both are used for b-cell localization to b/t border

globals [                          ;; Both globals below are used to measure the time it takes for an infection to clear

  naive-b-cell-spawn-timer
  llpc-lifespan-timer
  num-bacteria-escaped-to-blood
  num-naive-b-cells-in-blood
  num-llpcs-in-blood
  is-gc-seeded

  days-passed

  counter
]


;Called when the "setup" button is clicked. Should be the first action by the user.
to setup
  clear-all

  set days-passed 0

  if recording [
    carefully [ vid:start-recorder ] [ user-message error-message ]
  ]


  ;; Sets up the world structure (lymph node follicle + surrounding paracortex)
  ask patch 0 0 [ask patches in-radius 200  [set patch-type 1 set pcolor gray ]] ; outer zone
  ;ask patch 50 0 [ set pcolor gray set patch-type 1 ]
  ask patch 0 0 [ask patches in-radius 49  [set patch-type 0 set pcolor black]] ; inner zone
  ask patches with [ (pxcor = -50 or pxcor = -49) and abs(pycor) < 5 ] [ set patch-type 2 set pcolor red ]

  ask patches [ set s1p-level 0 set cxcl13-level 0 set ccl19-level 0 set ebi2-level 0 ]

  create-fdcs 60             ;; creates arbitrary # of FDCs for now
  ask fdcs [ set shape "square" set color brown setxy (-20 + random 40) (-20 + random 40) ] ;;non-possible antigen value during spawn to represent FDC with no antigen bound

  create-tfh-cells 20
  ask tfh-cells [ move-to one-of patches with [patch-type = 1] set time-alive -1000 set shape "square" set color cyan  set cxcr5-level 10 set ccr7-level 6 set ebi2r-level 6 set bcell-binding-status false]


  ;; Initialize global variables and counters
  set counter 0

  set num-naive-b-cells-in-blood 1000
  set is-gc-seeded false



  reset-ticks
end



;Called every tick
to go

  if recording [
    vid:record-interface

  ]

  set days-passed counter / 48

  spawn-b-cell
  spawn-th0-cell


  ask patches with [ patch-type = 1 ] [
   set ccl19-level ccl19-level + 2
   set ebi2-level ebi2-level + 2
  ]
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
  ask tfh-cells [ tfh-cell-function ]
  ask th0-cells [th0-cell-function ]
  ask th1-cells [th1-cell-function ]
  ask th2-cells [th2-cell-function ]
  ask bacteria [ bacteria-function ]
  ask breg-cells [ breg-function ]

  ask antibodies [antibodies-function]

  update-chemokine-gradient

  insert-cytokines

  incoming-sepsis


;  ifelse counter > 100
;  [set counter 0][set counter counter + 1]
  set counter counter + 1

  tick
end

to incoming-sepsis
  ask patches [set tnf-a tnf-a + ((count bacteria) / 1000)]
  ask patches [set il6 il6 + ((count bacteria) / 450)]

end

;to end-sepsis
;  set incoming-tnf-a 0
;  set incoming-il6 0
;end

to spawn-b-cell

  if counter mod 10 = 0 and num-naive-b-cells-in-blood != 0[
    create-naive-b-cells 1 [ set shape "circle" set color white set size 1 setxy 49 0
      set time-alive 0
      set bcr random 30
      set isotype "md"              ;; isotype of "md" is IgM/IgD coexpresion. "d" is IgD, "m" is IgM, "a" is IgA, "g" is IgG, "e" is IgE
      set s1pr1-level 0
      set s1pr2-level 0
      set cxcr5-level 16
      set ccr7-level 0
      set ebi2r-level 0
      set cd21-level 100

      ;set num-naive-b-cells-in-blood num-naive-b-cells-in-blood - 1
      set in-blood false
    ]
  ]
end

to spawn-th0-cell

  if counter mod 20 = 0 [
    create-th0-cells 1 [ set shape "square" set color yellow
      move-to one-of patches with [patch-type = 1]
      set time-alive 0
      set s1pr1-level 0
      set s1pr2-level 0
      set cxcr5-level 0
      set ccr7-level 6
      set ebi2r-level 6
      ;set num-naive-b-cells-in-blood num-naive-b-cells-in-blood - 1
      set in-blood false
    ]
  ]
end


to fdc-function
  set cxcl13-level cxcl13-level + 2
  set il4 il4 + 1
  set il6 il6 + 2
  set il15 il15 + 2
  set il12 il12 + 2

  if presented-antigen != 0 [
    set time-presenting time-presenting + 1
  ]
  if time-presenting = 300 [
    set presented-antigen 0
    set color brown
    set presented-antigen-type 0
    set time-presenting 0
  ]
end

to antibodies-function
   set time-alive time-alive + 1
  if time-alive > 1200 [
    die
  ]
end

to naive-b-cell-function

  set cd21-level 20 - ((il6 - il10) * 10)

  if in-blood = false [
    if patch-type = 2 [ ;; naive b cell exits LN
      set in-blood true
      ;hide-turtle
      die
    ]


if cd21-level > 10[
    ;; First, checks if naive b-cell is in contact with an APC presenting an antigen
    ;; Selects a single bacteria/antigen at the naive-b-cells current location
    let apc one-of fdcs-here
    let antigen one-of bacteria-here
    if (apc != nobody and [presented-antigen] of apc != 0) or antigen != nobody[
      set breed activated-b-cells
      set pro-breg 0
      set shape "circle"
      set size 1
      set color yellow
      set csr-bool false
      set time-alive 0
      ifelse antigen != nobody [
        let rTI random [num-TI-ag] of antigen
        let rTD random [num-TD-ag] of antigen
        ifelse rTI > rTD [
          set response-type 1
        ][
          set response-type 2   ; 2 is TD
          set ccr7-level 12
          set ebi2r-level 12
        ]

        ask antigen [ die ]


      ][
        if apc != nobody [
          set response-type [presented-antigen-type] of apc
        ]
      ]
    ]
    ]
    check-breg-status

    chemotaxis
    move

    if time-alive > 300 [
      set s1pr1-level s1pr1-level + 0.5 ;; this slowly increases the # of s1p receptors (s1pr) in the naive b cell when the b-cell is old enough
    ]
  ]

  check-tnf-status

  set time-alive time-alive + 1
  if time-alive > 1000 [
    die
  ]

end

to check-breg-status
  ifelse pro-breg > 250 [
    set breed breg-cells
    set size 1
    set shape "circle"
    set color violet
    set s1pr1-level 0
    set time-alive 0
  ][
    set pro-breg pro-breg + il6 + il21
  ]



end

to check-tnf-status
  set tnf-a-stimulation 100 * tnf-a

  if tnf-a-stimulation > 60 [
    ;print "apoptose"
    let x random 900
    if x = 0 [
      die
    ]
  ]
end

to breg-function
  set il10 il10 + 5
  set tgf-b tgf-b + 1
  chemotaxis
  move

  check-tnf-status

  set time-alive time-alive + 1
  if time-alive > 300 [
    die
  ]
end

to activated-b-cell-function
  if in-blood = false [
    if patch-type = 2 [
      set in-blood true
      hide-turtle
    ]

    isotype-switch

    ifelse response-type = 2 [
      td-response
    ][
      ifelse response-type = 1 [
        ti-response
      ][
        activated-mem-response
      ]
    ]

    check-breg-status

    chemotaxis
    move

  ]

  check-tnf-status

  set time-alive time-alive + 1
  if time-alive > 300 [
    die
  ]
end

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

to td-response
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
  ;; Activated B-cell for gc response upregulates CCR7 and EBI2R levels (capped out here at 12 for reasonably realistic localization behavior)
;  if ccr7-level < 12 [ set ccr7-level ccr7-level + 0.5 ]
;  if ebi2r-level < 12 [ set ebi2r-level ebi2r-level + 0.5 ]
end

to ti-response
  set time-alive time-alive
  set tnf-a tnf-a + 1
  if counter mod 20 = 0 [
    let proPC (il21 + il10 + if-a + if-g )
      let proMem (il21 + il4); * 100
    ;let rPC random proPC
    ;let rMem random proMem

    ifelse proPC > proMem [
      hatch-sl-plasma-cells 1 [ set time-alive 0 set color lime + 3 set shape "circle" set size 1 set s1pr1-level 0 set pro-breg 0]
    ][
      hatch-mem-b-cells 1 [set time-alive 0 set color white set shape "target" set s1pr1-level 10 set pro-breg 0 set cd21-level 100 set cxcr5-level 0]
    ]

  ]
end

to activated-mem-response
  set tnf-a tnf-a + 1
  if counter mod 100 = 0 [
    hatch-sl-plasma-cells 1 [ set time-alive 0 set color lime + 3 set shape "circle" set size 1 set s1pr1-level 0 set pro-breg 0]
  ]
end


to gc-b-cell-function

  if in-blood = false [
    if patch-type = 2 [
      set in-blood true
      hide-turtle
    ]

    check-breg-status

    set ebi2r-level 0
    set ccr7-level 0

    ifelse distance patch 0 0 > 15 [
      chemotaxis
      gc-move
    ][
      let proPC (il21 + il10 + if-a + if-g )
      let proMem (il21 + il4) * 2;* 100
;     let rPC random proPC
;     let rMem random proMem

      set level-of-activation il2 + il4 + il10 + il15 + il21 - if-g - if-a
      ;if round (time-alive mod (20 / level-of-activation)) = 0 [
      if time-alive mod 20 = 0 [
        ifelse proPC > proMem [
          hatch-ll-plasma-cells 1 [ set time-alive 0 set color lime set shape "circle" set size 1 set s1pr1-level 40 set pro-breg 0]
        ][
          hatch-mem-b-cells 1 [ set time-alive 0 set color white set shape "target" set s1pr1-level 10 set cxcr5-level 0 set pro-breg 0 set cd21-level 100]
        ]
      ]
    ]
  ]

  check-tnf-status

  set time-alive time-alive + 1
  if time-alive > 1000 [
    ask link-neighbors [ set bcell-binding-status false ]
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

  set level-of-activation il6
    ;if round (time-alive mod (10 / level-of-activation)) = 0 [
    if time-alive mod 50 = 0 [
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
      ;set num-llpcs-in-blood num-llpcs-in-blood + 1
      set in-blood true
      hide-turtle
    ]

    check-breg-status

    chemotaxis
    move
  ][
    set level-of-activation il6
    ;if round (time-alive mod (50 / level-of-activation)) = 0 [
    if time-alive mod 200 = 0 [
      hatch-antibodies 1 [ set time-alive 0 set antibody-type isotype set hidden? true  ]
    ]
  ]


  check-tnf-status

  set time-alive time-alive + 1
  if time-alive > 8000 + (il6 + il21) * 10 [
      die
  ]


end

to mem-b-cell-function
  set cd21-level 20 - ((il6 - il10) * 10)
  ifelse in-blood = false [
     if patch-type = 2 [
      ;set num-llpcs-in-blood num-llpcs-in-blood + 1
      set in-blood true
      hide-turtle
      set time-in-blood 0
    ]

;    let apc one-of fdcs-here
;    let antigen one-of bacteria-here
;    if (apc != nobody and [presented-antigen] of apc != 0) or antigen != nobody[
;      set breed activated-b-cells
;      set pro-breg 0
;      set shape "circle"
;      set size 1
;      set color yellow
;      set csr-bool false
;      set time-alive 0
;      set s1pr1-level 0
;      ifelse antigen != nobody [
;        set response-type 3  ;; 3 response type is memb reponse upon activation
;        ask antigen [ die ]
;
;      ][
;        if apc != nobody [
;          set response-type 3
;        ]
;      ]
;    ]

    if cd21-level > 10[
    let apc one-of fdcs-here
    let antigen one-of bacteria-here
    if (apc != nobody and [presented-antigen] of apc != 0) or antigen != nobody[
      set breed activated-b-cells
      set pro-breg 0
      set shape "circle"
      set size 1
      set color yellow
      set csr-bool false
      set time-alive 100
      ifelse antigen != nobody [
        let rTI random [num-TI-ag] of antigen
        let rTD random [num-TD-ag] of antigen
        ifelse rTI > rTD [
          set response-type 1
        ][
          set response-type 2   ; 2 is TD
          set ccr7-level 12
          set ebi2r-level 12
        ]

        ask antigen [ die ]


      ][
        if apc != nobody [
          set response-type [presented-antigen-type] of apc
        ]
      ]
    ]
    ]

    check-breg-status

    chemotaxis
    move
  ][
    set time-in-blood time-in-blood + 1
    if time-in-blood > 300 [
      if hidden? [
        set pro-breg 0
        set hidden? false
        set in-blood false
        setxy 49 0
        set s1pr1-level 10
        set cxcr5-level 10
      ]
    ]
  ]

  check-tnf-status



  set time-alive time-alive + 1
  if time-alive > 15000 [
    die
  ]

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

  set il21 il21 + 2
  set il4 il4 + 1
  set il2 il2 + 1
  set il10 il10 + 1

  set time-alive time-alive + 1
  if time-alive > 500
    [die]


end

to th1-cell-function
  chemotaxis
  move

  set if-g if-g + 1

  set time-alive time-alive + 1
  if time-alive > 500
    [die]

end

to th2-cell-function
  if distance patch 0 0 > 20 or bcell-binding-status = false [
    chemotaxis
    move
  ]

  set il4 il4 + 1
  set il10 il10 + 1

  set time-alive time-alive + 1
  if time-alive > 500 [
    die
  ]

end

to bacteria-function
  if patch-type = 2 [ ;; for bacteria, im having them recirculate through blood. when recirculating, they can either just go back into LN, or can be captured by FDC. random cahnce of either
    let x random 2

    ifelse x = 0 [
      setxy 49 0
    ][
      if any? fdcs with [presented-antigen = 0] [
        ask one-of fdcs with [presented-antigen = 0] [
          set time-presenting 0
          set presented-antigen bacteria-epitope-type
          set color red
          let rTI random number-of-TI-epitopes
          let rTD random number-of-TD-epitopes
          ifelse rTI > rTD [
            set presented-antigen-type 1   ;; 1 is TI
          ][
            set presented-antigen-type 2    ;; 2 is TD
          ]
        ]
      ]
    ]
  ]
  chemotaxis
  move


  set time-alive time-alive + 1
  ;if time-alive > 500
    ;[die]

end

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
    set il6 il6 * 0.9
    set il10 il10 * 0.9
    set il12 il12 * 0.9
    set il15 il15 * 0.9
    set il21 il21 * 0.9
    set if-g if-g * 0.9
    set if-a if-a * 0.9
    set tnf-a tnf-a * 0.9
    set tgf-b tgf-b * 0.9

    if patch-type = 0 [
      let total-cytokine-level il2 + il4 + il6 + il10 + il12 + il15 + il21 + tnf-a + tgf-b + if-a + if-g
      if cytokine != "none" [
        if cytokine = "tnf-a" [
          set pcolor scale-color green tnf-a 0.1 3  ;;used to visualize cxcl13 or ccl19 gradient
        ]
        if cytokine = "il6" [
          set pcolor scale-color green il6 0.1 3  ;;used to visualize cxcl13 or ccl19 gradient
        ]
        if cytokine = "il10" [
          set pcolor scale-color green il10 0.1 3  ;;used to visualize cxcl13 or ccl19 gradient
        ]
      ]
    ]


  ]
end




;This function is called when the user clicks the "inoculate" button in the interface. It adds bacteria into the system
to inoculate

  ;ask patches [ set tnf-a tnf-a + (number-of-bacteria / 10) ]
  ;set incoming-tnf-a (number-of-bacteria / 1000)
  ;set incoming-il6 (number-of-bacteria / 1000)
  ;ask patches [ set il6 il6 + (number-of-bacteria / 4) ]



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
    set num-TI-ag number-of-TI-epitopes
    set num-TD-ag number-of-TD-epitopes
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

to gc-move
  rt random 50
  lt random 50
  fd 0.5
end


to make-movie
  if vid:recorder-status = "inactive" [
    user-message "The recorder is inactive. There is nothing to save."
    stop
  ]
  ; prompt user for movie location
  user-message (word
    "Choose a name for your movie file (the "
    ".mp4 extension will be automatically added).")
  let path user-new-file
  if not is-string? path [ stop ]  ; stop if user canceled
  ; export the movie
  carefully [
    vid:save-recording path
    user-message (word "Exported movie to " path ".")
  ] [
    user-message error-message
  ]
end

to insert-cytokines
;  if mouse-down?     ;; reports true or false to indicate whether mouse button is down
;    [
;      ask patch mouse-xcor mouse-ycor [
;        ask patches in-radius 3 [
;          set il4 il-4
;          set il10 il-10
;          set il12 il-12
;          set il15 il-15
;          set il21 il-21
;          set if-a ifn-alpha
;          set if-g ifn-gamma
;          set tgf-b tgf-beta
;        ]
;      ]
;    ]
end
@#$#@#$#@
GRAPHICS-WINDOW
683
10
1133
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
14
16
106
54
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
119
16
217
55
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
548
23
634
56
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
334
11
512
44
number-of-bacteria
number-of-bacteria
0
500
52.0
1
1
NIL
HORIZONTAL

SLIDER
329
46
518
79
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
324
81
529
114
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
327
116
527
149
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
244
158
656
439
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
"LLPCs" 1.0 0 -14439633 true "" "plot count ll-plasma-cells with [in-blood = true]"
"SLPCs" 1.0 0 -8330359 true "" "plot count sl-plasma-cells"
"Mem B-Cells" 1.0 0 -7500403 true "" "plot count mem-b-cells"
"Total B Lymphocytes" 1.0 0 -2674135 true "" "plot (count mem-b-cells) + (count ll-plasma-cells with [in-blood = true]) + (count sl-plasma-cells )"

PLOT
17
179
192
306
Antibody Response Curve
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
16
316
197
440
IL-10 Production
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

CHOOSER
684
477
822
522
cytokine
cytokine
"none" "tnf-a" "il6" "il10"
0

BUTTON
125
72
248
105
NIL
make-movie
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
15
71
117
104
recording
recording
1
1
-1000

MONITOR
25
122
188
167
Number of Days Elapsed
days-passed
1
1
11

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
