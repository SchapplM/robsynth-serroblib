% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRPR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
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
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
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
	% StartTime: 2019-10-10 10:11:39
	% EndTime: 2019-10-10 10:11:39
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (59->20), mult. (192->41), div. (0->0), fcn. (208->8), ass. (0->25)
	t224 = sin(pkin(11));
	t226 = cos(pkin(11));
	t227 = cos(pkin(6));
	t230 = cos(qJ(2));
	t237 = qJD(2) * t230;
	t228 = sin(qJ(2));
	t238 = qJD(2) * t228;
	t212 = (t224 * t238 - t226 * t237) * t227;
	t234 = t230 * t224 + t228 * t226;
	t215 = t234 * t227;
	t217 = -t224 * t237 - t226 * t238;
	t218 = t228 * t224 - t230 * t226;
	t229 = sin(qJ(1));
	t231 = cos(qJ(1));
	t240 = -t229 * t212 - t231 * t217 + (t215 * t231 - t218 * t229) * qJD(1);
	t225 = sin(pkin(6));
	t239 = qJD(1) * t225;
	t233 = qJD(2) * t234;
	t216 = t218 * qJD(2);
	t232 = t231 * t212 - t229 * t217 + (t215 * t229 + t218 * t231) * qJD(1);
	t214 = t218 * t227;
	t213 = t227 * t233;
	t211 = -t231 * t213 + t229 * t216 + (t214 * t229 - t231 * t234) * qJD(1);
	t210 = t229 * t213 + t231 * t216 + (t214 * t231 + t229 * t234) * qJD(1);
	t1 = [t232, t210, 0, 0, 0, 0; -t240, t211, 0, 0, 0, 0; 0, -t225 * t233, 0, 0, 0, 0; -t211, t240, 0, 0, 0, 0; t210, t232, 0, 0, 0, 0; 0, t225 * t216, 0, 0, 0, 0; -t229 * t239, 0, 0, 0, 0, 0; t231 * t239, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:40
	% EndTime: 2019-10-10 10:11:40
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (191->44), mult. (594->89), div. (0->0), fcn. (666->10), ass. (0->46)
	t382 = sin(pkin(11));
	t384 = cos(pkin(11));
	t385 = cos(pkin(6));
	t390 = cos(qJ(2));
	t402 = qJD(2) * t390;
	t387 = sin(qJ(2));
	t403 = qJD(2) * t387;
	t367 = (t382 * t403 - t384 * t402) * t385;
	t397 = t390 * t382 + t387 * t384;
	t373 = t397 * t385;
	t375 = t387 * t382 - t390 * t384;
	t391 = cos(qJ(1));
	t388 = sin(qJ(1));
	t404 = qJD(1) * t388;
	t374 = -t382 * t402 - t384 * t403;
	t405 = t388 * t374;
	t359 = -t373 * t404 + t405 + (-qJD(1) * t375 - t367) * t391;
	t361 = t391 * t373 - t388 * t375;
	t386 = sin(qJ(4));
	t389 = cos(qJ(4));
	t383 = sin(pkin(6));
	t399 = t383 * t404;
	t406 = t383 * t391;
	t408 = (-t361 * t389 + t386 * t406) * qJD(4) - t359 * t386 + t389 * t399;
	t407 = t383 * t388;
	t401 = qJD(4) * t386;
	t400 = qJD(4) * t389;
	t398 = qJD(1) * t406;
	t372 = t375 * t385;
	t360 = -t391 * t372 - t388 * t397;
	t362 = t388 * t372 - t391 * t397;
	t363 = -t388 * t373 - t391 * t375;
	t370 = t375 * t383;
	t395 = qJD(2) * t397;
	t394 = t375 * qJD(2);
	t392 = -t359 * t389 + t400 * t406 + (qJD(4) * t361 - t399) * t386;
	t357 = -t361 * qJD(1) + t388 * t367 + t391 * t374;
	t371 = t397 * t383;
	t368 = t385 * t395;
	t366 = qJD(2) * t370;
	t365 = t383 * t395;
	t358 = t362 * qJD(1) - t391 * t368 + t388 * t394;
	t356 = t360 * qJD(1) - t388 * t368 - t391 * t394;
	t355 = t386 * t398 + t357 * t389 + (-t363 * t386 + t389 * t407) * qJD(4);
	t354 = t389 * t398 - t357 * t386 + (-t363 * t389 - t386 * t407) * qJD(4);
	t1 = [t392, -t356 * t389 - t362 * t401, 0, t354, 0, 0; t355, t358 * t389 - t360 * t401, 0, t408, 0, 0; 0, -t365 * t389 + t370 * t401, 0, t366 * t386 + (-t371 * t389 - t385 * t386) * qJD(4), 0, 0; -t408, t356 * t386 - t362 * t400, 0, -t355, 0, 0; t354, -t358 * t386 - t360 * t400, 0, t392, 0, 0; 0, t365 * t386 + t370 * t400, 0, t366 * t389 + (t371 * t386 - t385 * t389) * qJD(4), 0, 0; t358, t357, 0, 0, 0, 0; t356, t363 * qJD(1) - t391 * t367 + t405, 0, 0, 0, 0; 0, -t366, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:41
	% EndTime: 2019-10-10 10:11:42
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (314->55), mult. (970->115), div. (0->0), fcn. (1082->12), ass. (0->53)
	t483 = sin(pkin(11));
	t486 = cos(pkin(11));
	t487 = cos(pkin(6));
	t492 = cos(qJ(2));
	t507 = qJD(2) * t492;
	t489 = sin(qJ(2));
	t508 = qJD(2) * t489;
	t465 = (t483 * t508 - t486 * t507) * t487;
	t502 = t492 * t483 + t489 * t486;
	t471 = t502 * t487;
	t473 = t489 * t483 - t492 * t486;
	t493 = cos(qJ(1));
	t490 = sin(qJ(1));
	t509 = qJD(1) * t490;
	t472 = -t483 * t507 - t486 * t508;
	t510 = t490 * t472;
	t456 = -t471 * t509 + t510 + (-qJD(1) * t473 - t465) * t493;
	t459 = t493 * t471 - t490 * t473;
	t488 = sin(qJ(4));
	t491 = cos(qJ(4));
	t484 = sin(pkin(6));
	t504 = t484 * t509;
	t511 = t484 * t493;
	t450 = (t459 * t488 + t491 * t511) * qJD(4) - t456 * t491 - t488 * t504;
	t512 = t484 * t490;
	t506 = qJD(4) * t488;
	t505 = qJD(4) * t491;
	t503 = qJD(1) * t511;
	t470 = t473 * t487;
	t458 = -t493 * t470 - t490 * t502;
	t460 = t490 * t470 - t493 * t502;
	t461 = -t490 * t471 - t493 * t473;
	t497 = qJD(2) * t502;
	t466 = t487 * t497;
	t496 = t473 * qJD(2);
	t452 = t458 * qJD(1) - t490 * t466 - t493 * t496;
	t500 = t452 * t491 + t460 * t506;
	t455 = t460 * qJD(1) - t493 * t466 + t490 * t496;
	t499 = -t455 * t491 + t458 * t506;
	t463 = t484 * t497;
	t468 = t473 * t484;
	t498 = t463 * t491 - t468 * t506;
	t449 = -t456 * t488 - t459 * t505 + t491 * t504 + t506 * t511;
	t494 = -t459 * qJD(1) + t490 * t465 + t493 * t472;
	t485 = cos(pkin(12));
	t482 = sin(pkin(12));
	t469 = t502 * t484;
	t464 = qJD(2) * t468;
	t457 = t464 * t488 + (-t469 * t491 - t487 * t488) * qJD(4);
	t454 = t461 * qJD(1) - t493 * t465 + t510;
	t448 = t488 * t503 + t494 * t491 + (-t461 * t488 + t491 * t512) * qJD(4);
	t447 = t494 * t488 - t491 * t503 + (t461 * t491 + t488 * t512) * qJD(4);
	t1 = [t450 * t485 + t455 * t482, t482 * t494 - t500 * t485, 0, -t447 * t485, 0, 0; t448 * t485 + t452 * t482, t454 * t482 - t499 * t485, 0, t449 * t485, 0, 0; 0, -t464 * t482 - t498 * t485, 0, t457 * t485, 0, 0; -t450 * t482 + t455 * t485, t500 * t482 + t485 * t494, 0, t447 * t482, 0, 0; -t448 * t482 + t452 * t485, t454 * t485 + t499 * t482, 0, -t449 * t482, 0, 0; 0, -t464 * t485 + t498 * t482, 0, -t457 * t482, 0, 0; t449, -t452 * t488 + t460 * t505, 0, t448, 0, 0; t447, t455 * t488 + t458 * t505, 0, -t450, 0, 0; 0, -t463 * t488 - t468 * t505, 0, -t464 * t491 + (-t469 * t488 + t487 * t491) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:42
	% EndTime: 2019-10-10 10:11:43
	% DurationCPUTime: 0.63s
	% Computational Cost: add. (637->83), mult. (1658->154), div. (0->0), fcn. (1910->12), ass. (0->70)
	t583 = sin(pkin(11));
	t585 = cos(pkin(11));
	t586 = cos(pkin(6));
	t591 = cos(qJ(2));
	t614 = qJD(2) * t591;
	t588 = sin(qJ(2));
	t615 = qJD(2) * t588;
	t562 = (t583 * t615 - t585 * t614) * t586;
	t600 = t591 * t583 + t588 * t585;
	t568 = t600 * t586;
	t570 = t588 * t583 - t591 * t585;
	t592 = cos(qJ(1));
	t589 = sin(qJ(1));
	t616 = qJD(1) * t589;
	t569 = -t583 * t614 - t585 * t615;
	t618 = t589 * t569;
	t541 = -t568 * t616 + t618 + (-qJD(1) * t570 - t562) * t592;
	t587 = sin(qJ(4));
	t590 = cos(qJ(4));
	t602 = t592 * t568 - t589 * t570;
	t584 = sin(pkin(6));
	t620 = t584 * t592;
	t599 = t587 * t602 + t590 * t620;
	t607 = t584 * t616;
	t535 = t599 * qJD(4) - t541 * t590 - t587 * t607;
	t563 = qJD(2) * t568;
	t567 = t570 * t586;
	t598 = t570 * qJD(2);
	t540 = t567 * t616 + (-qJD(1) * t600 - t563) * t592 + t589 * t598;
	t608 = t587 * t620;
	t546 = -t590 * t602 + t608;
	t550 = -t592 * t567 - t589 * t600;
	t582 = pkin(12) + qJ(6);
	t580 = sin(t582);
	t581 = cos(t582);
	t629 = t535 * t581 + t540 * t580 + (-t546 * t580 + t550 * t581) * qJD(6);
	t628 = (t546 * t581 + t550 * t580) * qJD(6) + t535 * t580 - t540 * t581;
	t621 = t584 * t589;
	t613 = qJD(4) * t587;
	t612 = qJD(4) * t590;
	t611 = qJD(6) * t580;
	t610 = qJD(6) * t581;
	t609 = qJD(6) * t590;
	t606 = qJD(1) * t620;
	t553 = t589 * t567 - t592 * t600;
	t593 = -t602 * qJD(1) + t589 * t562 + t592 * t569;
	t605 = -t553 * t609 + t593;
	t601 = -t589 * t568 - t592 * t570;
	t604 = t601 * qJD(1) - t550 * t609 - t592 * t562 + t618;
	t561 = t584 * t598;
	t565 = t570 * t584;
	t603 = t565 * t609 - t561;
	t566 = t600 * t584;
	t556 = t566 * t590 + t586 * t587;
	t555 = -t566 * t587 + t586 * t590;
	t547 = -t587 * t601 + t590 * t621;
	t548 = t587 * t621 + t590 * t601;
	t533 = qJD(4) * t608 - t541 * t587 + t590 * t607 - t602 * t612;
	t537 = t550 * qJD(1) - t589 * t563 - t592 * t598;
	t596 = -qJD(6) * t601 + t537 * t590 + t553 * t613;
	t595 = -qJD(6) * t602 - t540 * t590 + t550 * t613;
	t560 = qJD(2) * t566;
	t594 = qJD(6) * t566 - t560 * t590 + t565 * t613;
	t543 = t555 * qJD(4) - t561 * t590;
	t542 = -t556 * qJD(4) + t561 * t587;
	t532 = t547 * qJD(4) + t587 * t606 + t590 * t593;
	t531 = t548 * qJD(4) + t587 * t593 - t590 * t606;
	t530 = t532 * t581 + t537 * t580 + (-t548 * t580 - t553 * t581) * qJD(6);
	t529 = -t532 * t580 + t537 * t581 + (-t548 * t581 + t553 * t580) * qJD(6);
	t1 = [t629, t605 * t580 - t596 * t581, 0, -t531 * t581 - t547 * t611, 0, t529; t530, t604 * t580 - t595 * t581, 0, t533 * t581 + t599 * t611, 0, t628; 0, t603 * t580 + t594 * t581, 0, t542 * t581 - t555 * t611, 0, -t543 * t580 + t560 * t581 + (-t556 * t581 - t565 * t580) * qJD(6); -t628, t596 * t580 + t605 * t581, 0, t531 * t580 - t547 * t610, 0, -t530; t529, t595 * t580 + t604 * t581, 0, -t533 * t580 + t599 * t610, 0, t629; 0, -t594 * t580 + t603 * t581, 0, -t542 * t580 - t555 * t610, 0, -t543 * t581 - t560 * t580 + (t556 * t580 - t565 * t581) * qJD(6); t533, -t537 * t587 + t553 * t612, 0, t532, 0, 0; t531, t540 * t587 + t550 * t612, 0, -t535, 0, 0; 0, -t560 * t587 - t565 * t612, 0, t543, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end