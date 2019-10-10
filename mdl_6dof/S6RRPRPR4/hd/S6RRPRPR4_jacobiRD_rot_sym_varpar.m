% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR4
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
% Datum: 2019-10-10 10:09
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRPR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:41
	% EndTime: 2019-10-10 10:09:41
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
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:42
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
	% StartTime: 2019-10-10 10:09:42
	% EndTime: 2019-10-10 10:09:42
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
	% StartTime: 2019-10-10 10:09:44
	% EndTime: 2019-10-10 10:09:44
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
	% StartTime: 2019-10-10 10:09:44
	% EndTime: 2019-10-10 10:09:44
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (241->45), mult. (594->89), div. (0->0), fcn. (666->10), ass. (0->47)
	t397 = sin(pkin(11));
	t399 = cos(pkin(11));
	t400 = cos(pkin(6));
	t403 = cos(qJ(2));
	t415 = qJD(2) * t403;
	t401 = sin(qJ(2));
	t416 = qJD(2) * t401;
	t379 = (t397 * t416 - t399 * t415) * t400;
	t410 = t403 * t397 + t401 * t399;
	t385 = t410 * t400;
	t387 = t401 * t397 - t403 * t399;
	t404 = cos(qJ(1));
	t402 = sin(qJ(1));
	t417 = qJD(1) * t402;
	t386 = -t397 * t415 - t399 * t416;
	t418 = t402 * t386;
	t371 = -t385 * t417 + t418 + (-qJD(1) * t387 - t379) * t404;
	t373 = t404 * t385 - t402 * t387;
	t396 = qJ(4) + pkin(12);
	t394 = sin(t396);
	t395 = cos(t396);
	t398 = sin(pkin(6));
	t412 = t398 * t417;
	t419 = t398 * t404;
	t421 = (-t373 * t395 + t394 * t419) * qJD(4) - t371 * t394 + t395 * t412;
	t420 = t398 * t402;
	t414 = qJD(4) * t394;
	t413 = qJD(4) * t395;
	t411 = qJD(1) * t419;
	t384 = t387 * t400;
	t372 = -t404 * t384 - t402 * t410;
	t374 = t402 * t384 - t404 * t410;
	t375 = -t402 * t385 - t404 * t387;
	t382 = t387 * t398;
	t408 = qJD(2) * t410;
	t407 = t387 * qJD(2);
	t405 = -t371 * t395 + t413 * t419 + (qJD(4) * t373 - t412) * t394;
	t369 = -t373 * qJD(1) + t402 * t379 + t404 * t386;
	t383 = t410 * t398;
	t380 = t400 * t408;
	t378 = qJD(2) * t382;
	t377 = t398 * t408;
	t370 = t374 * qJD(1) - t404 * t380 + t402 * t407;
	t368 = t372 * qJD(1) - t402 * t380 - t404 * t407;
	t367 = t394 * t411 + t369 * t395 + (-t375 * t394 + t395 * t420) * qJD(4);
	t366 = t395 * t411 - t369 * t394 + (-t375 * t395 - t394 * t420) * qJD(4);
	t1 = [t405, -t368 * t395 - t374 * t414, 0, t366, 0, 0; t367, t370 * t395 - t372 * t414, 0, t421, 0, 0; 0, -t377 * t395 + t382 * t414, 0, t378 * t394 + (-t383 * t395 - t394 * t400) * qJD(4), 0, 0; -t421, t368 * t394 - t374 * t413, 0, -t367, 0, 0; t366, -t370 * t394 - t372 * t413, 0, t405, 0, 0; 0, t377 * t394 + t382 * t413, 0, t378 * t395 + (t383 * t394 - t395 * t400) * qJD(4), 0, 0; t370, t369, 0, 0, 0, 0; t368, t375 * qJD(1) - t404 * t379 + t418, 0, 0, 0, 0; 0, -t378, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:09:46
	% EndTime: 2019-10-10 10:09:47
	% DurationCPUTime: 0.65s
	% Computational Cost: add. (690->83), mult. (1658->154), div. (0->0), fcn. (1910->12), ass. (0->70)
	t590 = sin(pkin(11));
	t592 = cos(pkin(11));
	t593 = cos(pkin(6));
	t598 = cos(qJ(2));
	t621 = qJD(2) * t598;
	t595 = sin(qJ(2));
	t622 = qJD(2) * t595;
	t569 = (t590 * t622 - t592 * t621) * t593;
	t607 = t598 * t590 + t595 * t592;
	t575 = t607 * t593;
	t577 = t595 * t590 - t598 * t592;
	t599 = cos(qJ(1));
	t596 = sin(qJ(1));
	t623 = qJD(1) * t596;
	t576 = -t590 * t621 - t592 * t622;
	t625 = t596 * t576;
	t548 = -t575 * t623 + t625 + (-qJD(1) * t577 - t569) * t599;
	t589 = qJ(4) + pkin(12);
	t587 = sin(t589);
	t588 = cos(t589);
	t609 = t599 * t575 - t596 * t577;
	t591 = sin(pkin(6));
	t627 = t591 * t599;
	t606 = t587 * t609 + t588 * t627;
	t614 = t591 * t623;
	t542 = t606 * qJD(4) - t548 * t588 - t587 * t614;
	t570 = qJD(2) * t575;
	t574 = t577 * t593;
	t605 = t577 * qJD(2);
	t547 = t574 * t623 + (-qJD(1) * t607 - t570) * t599 + t596 * t605;
	t615 = t587 * t627;
	t553 = -t588 * t609 + t615;
	t557 = -t599 * t574 - t596 * t607;
	t594 = sin(qJ(6));
	t597 = cos(qJ(6));
	t636 = t542 * t597 + t547 * t594 + (-t553 * t594 + t557 * t597) * qJD(6);
	t635 = (t553 * t597 + t557 * t594) * qJD(6) + t542 * t594 - t547 * t597;
	t628 = t591 * t596;
	t620 = qJD(4) * t587;
	t619 = qJD(4) * t588;
	t618 = qJD(6) * t588;
	t617 = qJD(6) * t594;
	t616 = qJD(6) * t597;
	t613 = qJD(1) * t627;
	t560 = t596 * t574 - t599 * t607;
	t600 = -t609 * qJD(1) + t596 * t569 + t599 * t576;
	t612 = -t560 * t618 + t600;
	t608 = -t596 * t575 - t599 * t577;
	t611 = t608 * qJD(1) - t557 * t618 - t599 * t569 + t625;
	t568 = t591 * t605;
	t572 = t577 * t591;
	t610 = t572 * t618 - t568;
	t573 = t607 * t591;
	t563 = t573 * t588 + t593 * t587;
	t562 = -t573 * t587 + t593 * t588;
	t554 = -t587 * t608 + t588 * t628;
	t555 = t587 * t628 + t588 * t608;
	t540 = qJD(4) * t615 - t548 * t587 + t588 * t614 - t609 * t619;
	t544 = t557 * qJD(1) - t596 * t570 - t599 * t605;
	t603 = -qJD(6) * t608 + t544 * t588 + t560 * t620;
	t602 = -qJD(6) * t609 - t547 * t588 + t557 * t620;
	t567 = qJD(2) * t573;
	t601 = qJD(6) * t573 - t567 * t588 + t572 * t620;
	t550 = t562 * qJD(4) - t568 * t588;
	t549 = -t563 * qJD(4) + t568 * t587;
	t539 = t554 * qJD(4) + t587 * t613 + t588 * t600;
	t538 = t555 * qJD(4) + t587 * t600 - t588 * t613;
	t537 = t539 * t597 + t544 * t594 + (-t555 * t594 - t560 * t597) * qJD(6);
	t536 = -t539 * t594 + t544 * t597 + (-t555 * t597 + t560 * t594) * qJD(6);
	t1 = [t636, t612 * t594 - t603 * t597, 0, -t538 * t597 - t554 * t617, 0, t536; t537, t611 * t594 - t602 * t597, 0, t540 * t597 + t606 * t617, 0, t635; 0, t610 * t594 + t601 * t597, 0, t549 * t597 - t562 * t617, 0, -t550 * t594 + t567 * t597 + (-t563 * t597 - t572 * t594) * qJD(6); -t635, t603 * t594 + t612 * t597, 0, t538 * t594 - t554 * t616, 0, -t537; t536, t602 * t594 + t611 * t597, 0, -t540 * t594 + t606 * t616, 0, t636; 0, -t601 * t594 + t610 * t597, 0, -t549 * t594 - t562 * t616, 0, -t550 * t597 - t567 * t594 + (t563 * t594 - t572 * t597) * qJD(6); t540, -t544 * t587 + t560 * t619, 0, t539, 0, 0; t538, t547 * t587 + t557 * t619, 0, -t542, 0, 0; 0, -t567 * t587 - t572 * t619, 0, t550, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end