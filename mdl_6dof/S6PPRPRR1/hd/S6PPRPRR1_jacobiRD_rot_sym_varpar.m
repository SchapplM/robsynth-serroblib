% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PPRPRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:08
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PPRPRR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRPRR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobiRD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:44
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (20->14), mult. (80->39), div. (0->0), fcn. (88->10), ass. (0->20)
	t132 = sin(pkin(11));
	t138 = cos(pkin(6));
	t147 = t132 * t138;
	t133 = sin(pkin(7));
	t134 = sin(pkin(6));
	t146 = t133 * t134;
	t145 = t133 * t138;
	t135 = cos(pkin(12));
	t137 = cos(pkin(7));
	t144 = t135 * t137;
	t136 = cos(pkin(11));
	t143 = t136 * t138;
	t131 = sin(pkin(12));
	t142 = -(-t132 * t131 + t135 * t143) * t137 + t136 * t146;
	t141 = -(-t136 * t131 - t135 * t147) * t137 - t132 * t146;
	t140 = cos(qJ(3));
	t139 = sin(qJ(3));
	t130 = -t131 * t147 + t136 * t135;
	t128 = t131 * t143 + t132 * t135;
	t1 = [0, 0, (-t130 * t140 + t141 * t139) * qJD(3), 0, 0, 0; 0, 0, (-t128 * t140 + t142 * t139) * qJD(3), 0, 0, 0; 0, 0, (-t139 * t145 + (-t131 * t140 - t139 * t144) * t134) * qJD(3), 0, 0, 0; 0, 0, (t130 * t139 + t141 * t140) * qJD(3), 0, 0, 0; 0, 0, (t128 * t139 + t142 * t140) * qJD(3), 0, 0, 0; 0, 0, (-t140 * t145 + (t131 * t139 - t140 * t144) * t134) * qJD(3), 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:44
	% EndTime: 2019-10-09 21:08:44
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (38->18), mult. (146->42), div. (0->0), fcn. (154->12), ass. (0->27)
	t180 = sin(pkin(13));
	t185 = cos(pkin(13));
	t190 = sin(qJ(3));
	t191 = cos(qJ(3));
	t198 = (t180 * t191 + t185 * t190) * qJD(3);
	t182 = sin(pkin(11));
	t184 = sin(pkin(6));
	t197 = t182 * t184;
	t189 = cos(pkin(6));
	t196 = t182 * t189;
	t187 = cos(pkin(11));
	t195 = t184 * t187;
	t194 = t187 * t189;
	t178 = (t180 * t190 - t185 * t191) * qJD(3);
	t188 = cos(pkin(7));
	t186 = cos(pkin(12));
	t183 = sin(pkin(7));
	t181 = sin(pkin(12));
	t177 = -t181 * t196 + t187 * t186;
	t176 = -t187 * t181 - t186 * t196;
	t175 = t181 * t194 + t182 * t186;
	t174 = -t182 * t181 + t186 * t194;
	t173 = t188 * t198;
	t172 = t188 * t178;
	t171 = t183 * t198;
	t170 = t183 * t178;
	t1 = [0, 0, -t171 * t197 - t176 * t173 + t177 * t178, 0, 0, 0; 0, 0, t171 * t195 - t174 * t173 + t175 * t178, 0, 0, 0; 0, 0, -t189 * t171 + (-t173 * t186 + t178 * t181) * t184, 0, 0, 0; 0, 0, t170 * t197 + t176 * t172 + t177 * t198, 0, 0, 0; 0, 0, -t170 * t195 + t174 * t172 + t175 * t198, 0, 0, 0; 0, 0, t189 * t170 + (t172 * t186 + t181 * t198) * t184, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:45
	% EndTime: 2019-10-09 21:08:46
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (199->53), mult. (663->118), div. (0->0), fcn. (769->14), ass. (0->57)
	t410 = sin(pkin(13));
	t415 = cos(pkin(13));
	t423 = cos(qJ(3));
	t431 = qJD(3) * t423;
	t421 = sin(qJ(3));
	t432 = qJD(3) * t421;
	t438 = t410 * t432 - t415 * t431;
	t412 = sin(pkin(11));
	t414 = sin(pkin(6));
	t437 = t412 * t414;
	t419 = cos(pkin(6));
	t436 = t412 * t419;
	t417 = cos(pkin(11));
	t435 = t414 * t417;
	t418 = cos(pkin(7));
	t434 = t414 * t418;
	t433 = t417 * t419;
	t420 = sin(qJ(5));
	t430 = qJD(5) * t420;
	t422 = cos(qJ(5));
	t429 = qJD(5) * t422;
	t426 = t423 * t410 + t421 * t415;
	t425 = t421 * t410 - t423 * t415;
	t424 = qJD(3) * t426;
	t413 = sin(pkin(7));
	t389 = t438 * t413;
	t391 = t438 * t418;
	t411 = sin(pkin(12));
	t416 = cos(pkin(12));
	t398 = -t412 * t411 + t416 * t433;
	t399 = t411 * t433 + t412 * t416;
	t403 = -t410 * t431 - t415 * t432;
	t376 = t389 * t435 - t398 * t391 + t399 * t403;
	t400 = -t417 * t411 - t416 * t436;
	t401 = -t411 * t436 + t417 * t416;
	t378 = -t389 * t437 - t400 * t391 + t401 * t403;
	t384 = -t419 * t389 + (-t391 * t416 + t403 * t411) * t414;
	t402 = t425 * qJD(3);
	t397 = -t414 * t416 * t413 + t419 * t418;
	t396 = t426 * t418;
	t395 = t425 * t418;
	t394 = t426 * t413;
	t393 = t425 * t413;
	t392 = t418 * t424;
	t390 = t413 * t424;
	t388 = -t400 * t413 + t412 * t434;
	t387 = -t398 * t413 - t417 * t434;
	t386 = t419 * t394 + (t396 * t416 - t411 * t425) * t414;
	t385 = -t419 * t393 + (-t395 * t416 - t411 * t426) * t414;
	t383 = -t419 * t390 + (-t392 * t416 + t402 * t411) * t414;
	t382 = t394 * t437 + t400 * t396 - t401 * t425;
	t381 = -t393 * t437 - t400 * t395 - t401 * t426;
	t380 = -t394 * t435 + t398 * t396 - t399 * t425;
	t379 = t393 * t435 - t398 * t395 - t399 * t426;
	t377 = -t390 * t437 - t400 * t392 + t401 * t402;
	t375 = t390 * t435 - t398 * t392 + t399 * t402;
	t1 = [0, 0, t377 * t422 - t381 * t430, 0, -t378 * t420 + (-t382 * t422 - t388 * t420) * qJD(5), 0; 0, 0, t375 * t422 - t379 * t430, 0, -t376 * t420 + (-t380 * t422 - t387 * t420) * qJD(5), 0; 0, 0, t383 * t422 - t385 * t430, 0, -t384 * t420 + (-t386 * t422 - t397 * t420) * qJD(5), 0; 0, 0, -t377 * t420 - t381 * t429, 0, -t378 * t422 + (t382 * t420 - t388 * t422) * qJD(5), 0; 0, 0, -t375 * t420 - t379 * t429, 0, -t376 * t422 + (t380 * t420 - t387 * t422) * qJD(5), 0; 0, 0, -t383 * t420 - t385 * t429, 0, -t384 * t422 + (t386 * t420 - t397 * t422) * qJD(5), 0; 0, 0, t378, 0, 0, 0; 0, 0, t376, 0, 0, 0; 0, 0, t384, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:48
	% EndTime: 2019-10-09 21:08:49
	% DurationCPUTime: 0.60s
	% Computational Cost: add. (708->89), mult. (2229->181), div. (0->0), fcn. (2676->16), ass. (0->80)
	t584 = sin(pkin(13));
	t589 = cos(pkin(13));
	t599 = cos(qJ(3));
	t621 = qJD(3) * t599;
	t596 = sin(qJ(3));
	t622 = qJD(3) * t596;
	t628 = t584 * t622 - t589 * t621;
	t586 = sin(pkin(11));
	t588 = sin(pkin(6));
	t627 = t586 * t588;
	t593 = cos(pkin(6));
	t626 = t586 * t593;
	t591 = cos(pkin(11));
	t625 = t588 * t591;
	t592 = cos(pkin(7));
	t624 = t588 * t592;
	t623 = t591 * t593;
	t595 = sin(qJ(5));
	t620 = qJD(5) * t595;
	t598 = cos(qJ(5));
	t619 = qJD(5) * t598;
	t594 = sin(qJ(6));
	t618 = qJD(6) * t594;
	t597 = cos(qJ(6));
	t617 = qJD(6) * t597;
	t616 = qJD(6) * t598;
	t577 = t596 * t584 - t599 * t589;
	t587 = sin(pkin(7));
	t566 = t577 * t587;
	t568 = t577 * t592;
	t585 = sin(pkin(12));
	t590 = cos(pkin(12));
	t571 = -t586 * t585 + t590 * t623;
	t572 = t585 * t623 + t586 * t590;
	t610 = t599 * t584 + t596 * t589;
	t545 = t566 * t625 - t571 * t568 - t572 * t610;
	t562 = t628 * t587;
	t564 = t628 * t592;
	t576 = -t584 * t621 - t589 * t622;
	t608 = t562 * t625 - t571 * t564 + t572 * t576;
	t613 = -t545 * t616 + t608;
	t573 = -t591 * t585 - t590 * t626;
	t574 = -t585 * t626 + t591 * t590;
	t548 = -t566 * t627 - t573 * t568 - t574 * t610;
	t607 = -t562 * t627 - t573 * t564 + t574 * t576;
	t612 = -t548 * t616 + t607;
	t556 = -t593 * t566 + (-t568 * t590 - t585 * t610) * t588;
	t601 = -t593 * t562 + (-t564 * t590 + t576 * t585) * t588;
	t611 = -t556 * t616 + t601;
	t558 = -t571 * t587 - t591 * t624;
	t567 = t610 * t587;
	t569 = t610 * t592;
	t606 = -t567 * t625 + t571 * t569 - t572 * t577;
	t535 = t558 * t595 + t598 * t606;
	t534 = t558 * t598 - t595 * t606;
	t559 = -t573 * t587 + t586 * t624;
	t605 = t567 * t627 + t573 * t569 - t574 * t577;
	t537 = t559 * t595 + t598 * t605;
	t536 = t559 * t598 - t595 * t605;
	t570 = -t588 * t590 * t587 + t593 * t592;
	t600 = t593 * t567 + (t569 * t590 - t577 * t585) * t588;
	t551 = t570 * t595 + t598 * t600;
	t550 = t570 * t598 - t595 * t600;
	t609 = qJD(3) * t610;
	t563 = t587 * t609;
	t565 = t592 * t609;
	t575 = t577 * qJD(3);
	t538 = t563 * t625 - t571 * t565 + t572 * t575;
	t604 = -qJD(6) * t606 - t538 * t598 + t545 * t620;
	t541 = -t563 * t627 - t573 * t565 + t574 * t575;
	t603 = -qJD(6) * t605 - t541 * t598 + t548 * t620;
	t552 = -t593 * t563 + (-t565 * t590 + t575 * t585) * t588;
	t602 = -qJD(6) * t600 - t552 * t598 + t556 * t620;
	t533 = t550 * qJD(5) + t598 * t601;
	t532 = -t551 * qJD(5) - t595 * t601;
	t531 = t536 * qJD(5) + t598 * t607;
	t530 = -t537 * qJD(5) - t595 * t607;
	t529 = t534 * qJD(5) + t598 * t608;
	t528 = -t535 * qJD(5) - t595 * t608;
	t1 = [0, 0, t612 * t594 - t603 * t597, 0, t530 * t597 - t536 * t618, -t531 * t594 - t541 * t597 + (-t537 * t597 + t548 * t594) * qJD(6); 0, 0, t613 * t594 - t604 * t597, 0, t528 * t597 - t534 * t618, -t529 * t594 - t538 * t597 + (-t535 * t597 + t545 * t594) * qJD(6); 0, 0, t611 * t594 - t602 * t597, 0, t532 * t597 - t550 * t618, -t533 * t594 - t552 * t597 + (-t551 * t597 + t556 * t594) * qJD(6); 0, 0, t603 * t594 + t612 * t597, 0, -t530 * t594 - t536 * t617, -t531 * t597 + t541 * t594 + (t537 * t594 + t548 * t597) * qJD(6); 0, 0, t604 * t594 + t613 * t597, 0, -t528 * t594 - t534 * t617, -t529 * t597 + t538 * t594 + (t535 * t594 + t545 * t597) * qJD(6); 0, 0, t602 * t594 + t611 * t597, 0, -t532 * t594 - t550 * t617, -t533 * t597 + t552 * t594 + (t551 * t594 + t556 * t597) * qJD(6); 0, 0, t541 * t595 + t548 * t619, 0, t531, 0; 0, 0, t538 * t595 + t545 * t619, 0, t529, 0; 0, 0, t552 * t595 + t556 * t619, 0, t533, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end