% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:03
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRRR9_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR9_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (27->13), mult. (88->22), div. (0->0), fcn. (88->6), ass. (0->21)
	t158 = sin(qJ(2));
	t159 = sin(qJ(1));
	t171 = t158 * t159;
	t161 = cos(qJ(1));
	t170 = t158 * t161;
	t160 = cos(qJ(2));
	t169 = t159 * t160;
	t168 = t160 * t161;
	t156 = sin(pkin(6));
	t167 = qJD(1) * t156;
	t166 = qJD(2) * t156;
	t157 = cos(pkin(6));
	t165 = -t157 * t168 + t171;
	t164 = t157 * t169 + t170;
	t163 = t157 * t170 + t169;
	t162 = t157 * t171 - t168;
	t155 = t162 * qJD(1) + t165 * qJD(2);
	t154 = t164 * qJD(1) + t163 * qJD(2);
	t153 = t163 * qJD(1) + t164 * qJD(2);
	t152 = t165 * qJD(1) + t162 * qJD(2);
	t1 = [t155, t152, 0, 0, 0, 0; -t153, -t154, 0, 0, 0, 0; 0, -t158 * t166, 0, 0, 0, 0; t154, t153, 0, 0, 0, 0; t152, t155, 0, 0, 0, 0; 0, -t160 * t166, 0, 0, 0, 0; -t159 * t167, 0, 0, 0, 0, 0; t161 * t167, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:14
	% EndTime: 2019-10-10 11:03:14
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (43->18), mult. (148->39), div. (0->0), fcn. (148->8), ass. (0->25)
	t231 = sin(qJ(2));
	t232 = sin(qJ(1));
	t246 = t231 * t232;
	t234 = cos(qJ(1));
	t245 = t231 * t234;
	t233 = cos(qJ(2));
	t244 = t232 * t233;
	t243 = t233 * t234;
	t228 = sin(pkin(6));
	t242 = qJD(1) * t228;
	t241 = qJD(2) * t231;
	t230 = cos(pkin(6));
	t240 = t230 * t246;
	t239 = t232 * t242;
	t238 = t234 * t242;
	t237 = t228 * t241;
	t236 = -t230 * t244 - t245;
	t235 = -t230 * t245 - t244;
	t229 = cos(pkin(12));
	t227 = sin(pkin(12));
	t224 = -qJD(1) * t240 - t232 * t241 + (qJD(2) * t230 + qJD(1)) * t243;
	t223 = t236 * qJD(1) + t235 * qJD(2);
	t222 = t235 * qJD(1) + t236 * qJD(2);
	t221 = (t240 - t243) * qJD(2) + (-t230 * t243 + t246) * qJD(1);
	t1 = [-t224 * t229 - t227 * t239, t221 * t229, 0, 0, 0, 0; t222 * t229 + t227 * t238, t223 * t229, 0, 0, 0, 0; 0, -t229 * t237, 0, 0, 0, 0; t224 * t227 - t229 * t239, -t221 * t227, 0, 0, 0, 0; -t222 * t227 + t229 * t238, -t223 * t227, 0, 0, 0, 0; 0, t227 * t237, 0, 0, 0, 0; t223, t222, 0, 0, 0, 0; -t221, t224, 0, 0, 0, 0; 0, t228 * qJD(2) * t233, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:14
	% EndTime: 2019-10-10 11:03:15
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (144->36), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->38)
	t303 = cos(pkin(6));
	t305 = sin(qJ(1));
	t304 = sin(qJ(2));
	t324 = t305 * t304;
	t315 = t303 * t324;
	t319 = qJD(2) * t304;
	t306 = cos(qJ(2));
	t307 = cos(qJ(1));
	t321 = t307 * t306;
	t291 = -qJD(1) * t315 - t305 * t319 + (qJD(2) * t303 + qJD(1)) * t321;
	t322 = t307 * t304;
	t323 = t305 * t306;
	t293 = t303 * t322 + t323;
	t301 = pkin(12) + qJ(4);
	t299 = sin(t301);
	t300 = cos(t301);
	t302 = sin(pkin(6));
	t320 = qJD(1) * t302;
	t314 = t305 * t320;
	t325 = t302 * t307;
	t328 = (-t293 * t300 + t299 * t325) * qJD(4) - t291 * t299 + t300 * t314;
	t327 = t302 * t304;
	t326 = t302 * t305;
	t318 = qJD(4) * t299;
	t317 = qJD(4) * t300;
	t316 = qJD(4) * t306;
	t313 = t307 * t320;
	t312 = t302 * qJD(2) * t306;
	t292 = t303 * t321 - t324;
	t294 = -t303 * t323 - t322;
	t310 = t315 - t321;
	t308 = -t291 * t300 + t317 * t325 + (qJD(4) * t293 - t314) * t299;
	t290 = qJD(1) * t294 - qJD(2) * t293;
	t289 = -qJD(1) * t293 + qJD(2) * t294;
	t288 = -qJD(1) * t292 + qJD(2) * t310;
	t287 = t299 * t313 + t289 * t300 + (t299 * t310 + t300 * t326) * qJD(4);
	t286 = t300 * t313 - t289 * t299 + (-t299 * t326 + t300 * t310) * qJD(4);
	t1 = [t308, t288 * t300 - t294 * t318, 0, t286, 0, 0; t287, t290 * t300 - t292 * t318, 0, t328, 0, 0; 0, (-t299 * t316 - t300 * t319) * t302, 0, -t299 * t312 + (-t299 * t303 - t300 * t327) * qJD(4), 0, 0; -t328, -t288 * t299 - t294 * t317, 0, -t287, 0, 0; t286, -t290 * t299 - t292 * t317, 0, t308, 0, 0; 0, (t299 * t319 - t300 * t316) * t302, 0, -t300 * t312 + (t299 * t327 - t300 * t303) * qJD(4), 0, 0; t290, t289, 0, 0, 0, 0; -t288, t291, 0, 0, 0, 0; 0, t312, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:15
	% EndTime: 2019-10-10 11:03:15
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (308->33), mult. (414->59), div. (0->0), fcn. (430->8), ass. (0->43)
	t356 = cos(pkin(6));
	t358 = sin(qJ(1));
	t357 = sin(qJ(2));
	t378 = t358 * t357;
	t369 = t356 * t378;
	t372 = qJD(2) * t357;
	t359 = cos(qJ(2));
	t360 = cos(qJ(1));
	t374 = t360 * t359;
	t343 = -qJD(1) * t369 - t358 * t372 + (qJD(2) * t356 + qJD(1)) * t374;
	t354 = qJD(4) + qJD(5);
	t355 = sin(pkin(6));
	t379 = t354 * t355;
	t383 = t360 * t379 - t343;
	t364 = t369 - t374;
	t373 = qJD(1) * t355;
	t382 = t364 * t354 + t360 * t373;
	t353 = pkin(12) + qJ(4) + qJ(5);
	t351 = sin(t353);
	t352 = cos(t353);
	t375 = t360 * t357;
	t377 = t358 * t359;
	t347 = t356 * t375 + t377;
	t363 = -t347 * t354 + t358 * t373;
	t336 = t383 * t351 + t363 * t352;
	t381 = t351 * t354;
	t380 = t352 * t354;
	t376 = t359 * t354;
	t371 = t357 * t379;
	t367 = t355 * qJD(2) * t359;
	t348 = -t356 * t377 - t375;
	t341 = -t347 * qJD(1) + t348 * qJD(2);
	t366 = t358 * t379 + t341;
	t346 = t356 * t374 - t378;
	t362 = -t354 * t356 - t367;
	t337 = -t363 * t351 + t383 * t352;
	t342 = t348 * qJD(1) - t347 * qJD(2);
	t340 = -t346 * qJD(1) + t364 * qJD(2);
	t339 = t351 * t371 + t362 * t352;
	t338 = t362 * t351 - t352 * t371;
	t335 = t382 * t351 + t366 * t352;
	t334 = -t366 * t351 + t382 * t352;
	t1 = [t337, t340 * t352 - t348 * t381, 0, t334, t334, 0; t335, t342 * t352 - t346 * t381, 0, t336, t336, 0; 0, (-t351 * t376 - t352 * t372) * t355, 0, t338, t338, 0; -t336, -t340 * t351 - t348 * t380, 0, -t335, -t335, 0; t334, -t342 * t351 - t346 * t380, 0, t337, t337, 0; 0, (t351 * t372 - t352 * t376) * t355, 0, t339, t339, 0; t342, t341, 0, 0, 0, 0; -t340, t343, 0, 0, 0, 0; 0, t367, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:17
	% EndTime: 2019-10-10 11:03:18
	% DurationCPUTime: 0.67s
	% Computational Cost: add. (782->77), mult. (1102->143), div. (0->0), fcn. (1184->10), ass. (0->77)
	t552 = sin(qJ(1));
	t549 = cos(pkin(6));
	t566 = qJD(2) * t549 + qJD(1);
	t551 = sin(qJ(2));
	t585 = t552 * t551;
	t572 = t549 * t585;
	t580 = qJD(2) * t551;
	t554 = cos(qJ(2));
	t555 = cos(qJ(1));
	t582 = t555 * t554;
	t524 = -qJD(1) * t572 - t552 * t580 + t566 * t582;
	t583 = t555 * t551;
	t584 = t552 * t554;
	t534 = t549 * t583 + t584;
	t546 = pkin(12) + qJ(4) + qJ(5);
	t544 = sin(t546);
	t545 = cos(t546);
	t547 = qJD(4) + qJD(5);
	t548 = sin(pkin(6));
	t581 = qJD(1) * t548;
	t569 = t552 * t581;
	t586 = t548 * t555;
	t515 = (-t534 * t547 + t569) * t544 - (t547 * t586 - t524) * t545;
	t535 = t549 * t584 + t583;
	t523 = t535 * qJD(1) + t534 * qJD(2);
	t571 = t544 * t586;
	t527 = -t534 * t545 + t571;
	t570 = t549 * t582;
	t533 = -t570 + t585;
	t550 = sin(qJ(6));
	t553 = cos(qJ(6));
	t601 = -t515 * t553 + (-t527 * t550 - t533 * t553) * qJD(6) - t523 * t550;
	t600 = (t527 * t553 - t533 * t550) * qJD(6) - t515 * t550 + t523 * t553;
	t590 = t547 * t554;
	t597 = (qJD(2) * t545 - qJD(6)) * t551 + t544 * t590;
	t592 = t544 * t547;
	t591 = t545 * t547;
	t589 = t548 * t551;
	t588 = t548 * t552;
	t587 = t548 * t554;
	t579 = qJD(2) * t554;
	t578 = qJD(6) * t545;
	t577 = qJD(6) * t550;
	t576 = qJD(6) * t553;
	t574 = t544 * t589;
	t573 = t545 * t589;
	t568 = t555 * t581;
	t567 = t548 * t580;
	t522 = -t534 * qJD(1) - t535 * qJD(2);
	t564 = t547 * t588 + t522;
	t562 = t535 * t578 + t522;
	t561 = t533 * t578 + t524;
	t560 = (qJD(2) - t578) * t554;
	t558 = t547 * t549 + t548 * t579;
	t514 = -t524 * t544 - t534 * t591 + t545 * t569 + t547 * t571;
	t521 = -qJD(1) * t570 - t555 * t579 + t566 * t585;
	t536 = -t572 + t582;
	t557 = qJD(6) * t536 + t521 * t545 + t535 * t592;
	t556 = qJD(6) * t534 - t523 * t545 + t533 * t592;
	t531 = t549 * t544 + t573;
	t530 = t549 * t545 - t574;
	t529 = t536 * t545 + t544 * t588;
	t528 = -t536 * t544 + t545 * t588;
	t525 = -t534 * t544 - t545 * t586;
	t520 = t558 * t545 - t547 * t574;
	t519 = -t558 * t544 - t547 * t573;
	t518 = t519 * t553 - t530 * t577;
	t517 = -t519 * t550 - t530 * t576;
	t513 = t564 * t545 + (-t536 * t547 + t568) * t544;
	t512 = t536 * t591 + t564 * t544 - t545 * t568;
	t511 = t514 * t553 - t525 * t577;
	t510 = -t514 * t550 - t525 * t576;
	t509 = -t512 * t553 - t528 * t577;
	t508 = t512 * t550 - t528 * t576;
	t507 = t513 * t553 - t521 * t550 + (-t529 * t550 + t535 * t553) * qJD(6);
	t506 = -t513 * t550 - t521 * t553 + (-t529 * t553 - t535 * t550) * qJD(6);
	t1 = [t601, t562 * t550 + t557 * t553, 0, t509, t509, t506; t507, t561 * t550 + t556 * t553, 0, t511, t511, t600; 0, (t550 * t560 - t597 * t553) * t548, 0, t518, t518, t553 * t567 - t520 * t550 + (-t531 * t553 + t550 * t587) * qJD(6); -t600, -t557 * t550 + t562 * t553, 0, t508, t508, -t507; t506, -t556 * t550 + t561 * t553, 0, t510, t510, t601; 0, (t597 * t550 + t553 * t560) * t548, 0, t517, t517, -t550 * t567 - t520 * t553 + (t531 * t550 + t553 * t587) * qJD(6); t514, t521 * t544 - t535 * t591, 0, t513, t513, 0; t512, -t523 * t544 - t533 * t591, 0, t515, t515, 0; 0, (-t544 * t580 + t545 * t590) * t548, 0, t520, t520, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end