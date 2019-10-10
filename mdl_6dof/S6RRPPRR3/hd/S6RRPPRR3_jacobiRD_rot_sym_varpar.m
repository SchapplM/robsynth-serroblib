% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:39
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPPRR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:19
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
	% StartTime: 2019-10-10 09:39:19
	% EndTime: 2019-10-10 09:39:20
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
	% StartTime: 2019-10-10 09:39:20
	% EndTime: 2019-10-10 09:39:20
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (94->25), mult. (304->55), div. (0->0), fcn. (328->10), ass. (0->30)
	t317 = sin(pkin(11));
	t320 = cos(pkin(11));
	t322 = sin(qJ(2));
	t324 = cos(qJ(2));
	t310 = t322 * t317 - t324 * t320;
	t326 = t310 * qJD(2);
	t321 = cos(pkin(6));
	t333 = qJD(2) * t324;
	t334 = qJD(2) * t322;
	t304 = (t317 * t334 - t320 * t333) * t321;
	t328 = t324 * t317 + t322 * t320;
	t308 = t328 * t321;
	t309 = -t317 * t333 - t320 * t334;
	t323 = sin(qJ(1));
	t325 = cos(qJ(1));
	t302 = t325 * t304 - t323 * t309 + (t308 * t323 + t310 * t325) * qJD(1);
	t318 = sin(pkin(6));
	t335 = qJD(1) * t318;
	t332 = t323 * t335;
	t331 = t325 * t335;
	t327 = qJD(2) * t328;
	t300 = t323 * t304 + t325 * t309 + (-t308 * t325 + t310 * t323) * qJD(1);
	t319 = cos(pkin(12));
	t316 = sin(pkin(12));
	t307 = t310 * t321;
	t305 = t321 * t327;
	t303 = t318 * t327;
	t301 = -t325 * t305 + t323 * t326 + (t323 * t307 - t325 * t328) * qJD(1);
	t299 = -t323 * t305 - t325 * t326 + (-t307 * t325 - t323 * t328) * qJD(1);
	t1 = [t302 * t319 - t316 * t332, -t299 * t319, 0, 0, 0, 0; t300 * t319 + t316 * t331, t301 * t319, 0, 0, 0, 0; 0, -t303 * t319, 0, 0, 0, 0; -t302 * t316 - t319 * t332, t299 * t316, 0, 0, 0, 0; -t300 * t316 + t319 * t331, -t301 * t316, 0, 0, 0, 0; 0, t303 * t316, 0, 0, 0, 0; t301, t300, 0, 0, 0, 0; t299, -t302, 0, 0, 0, 0; 0, -t318 * t326, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:21
	% EndTime: 2019-10-10 09:39:21
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (241->45), mult. (594->89), div. (0->0), fcn. (666->10), ass. (0->47)
	t396 = sin(pkin(11));
	t398 = cos(pkin(11));
	t399 = cos(pkin(6));
	t402 = cos(qJ(2));
	t414 = qJD(2) * t402;
	t400 = sin(qJ(2));
	t415 = qJD(2) * t400;
	t378 = (t396 * t415 - t398 * t414) * t399;
	t409 = t402 * t396 + t400 * t398;
	t384 = t409 * t399;
	t386 = t400 * t396 - t402 * t398;
	t403 = cos(qJ(1));
	t401 = sin(qJ(1));
	t416 = qJD(1) * t401;
	t385 = -t396 * t414 - t398 * t415;
	t417 = t401 * t385;
	t370 = -t384 * t416 + t417 + (-qJD(1) * t386 - t378) * t403;
	t372 = t403 * t384 - t401 * t386;
	t395 = pkin(12) + qJ(5);
	t393 = sin(t395);
	t394 = cos(t395);
	t397 = sin(pkin(6));
	t411 = t397 * t416;
	t418 = t397 * t403;
	t420 = (-t372 * t394 + t393 * t418) * qJD(5) - t370 * t393 + t394 * t411;
	t419 = t397 * t401;
	t413 = qJD(5) * t393;
	t412 = qJD(5) * t394;
	t410 = qJD(1) * t418;
	t383 = t386 * t399;
	t371 = -t403 * t383 - t401 * t409;
	t373 = t401 * t383 - t403 * t409;
	t374 = -t401 * t384 - t403 * t386;
	t381 = t386 * t397;
	t407 = qJD(2) * t409;
	t406 = t386 * qJD(2);
	t404 = -t370 * t394 + t412 * t418 + (qJD(5) * t372 - t411) * t393;
	t368 = -t372 * qJD(1) + t401 * t378 + t403 * t385;
	t382 = t409 * t397;
	t379 = t399 * t407;
	t377 = qJD(2) * t381;
	t376 = t397 * t407;
	t369 = t373 * qJD(1) - t403 * t379 + t401 * t406;
	t367 = t371 * qJD(1) - t401 * t379 - t403 * t406;
	t366 = t393 * t410 + t368 * t394 + (-t374 * t393 + t394 * t419) * qJD(5);
	t365 = t394 * t410 - t368 * t393 + (-t374 * t394 - t393 * t419) * qJD(5);
	t1 = [t404, -t367 * t394 - t373 * t413, 0, 0, t365, 0; t366, t369 * t394 - t371 * t413, 0, 0, t420, 0; 0, -t376 * t394 + t381 * t413, 0, 0, t377 * t393 + (-t382 * t394 - t393 * t399) * qJD(5), 0; -t420, t367 * t393 - t373 * t412, 0, 0, -t366, 0; t365, -t369 * t393 - t371 * t412, 0, 0, t404, 0; 0, t376 * t393 + t381 * t412, 0, 0, t377 * t394 + (t382 * t393 - t394 * t399) * qJD(5), 0; t369, t368, 0, 0, 0, 0; t367, t374 * qJD(1) - t403 * t378 + t417, 0, 0, 0, 0; 0, -t377, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:23
	% EndTime: 2019-10-10 09:39:24
	% DurationCPUTime: 0.63s
	% Computational Cost: add. (690->83), mult. (1658->154), div. (0->0), fcn. (1910->12), ass. (0->70)
	t588 = sin(pkin(11));
	t590 = cos(pkin(11));
	t591 = cos(pkin(6));
	t596 = cos(qJ(2));
	t619 = qJD(2) * t596;
	t593 = sin(qJ(2));
	t620 = qJD(2) * t593;
	t567 = (t588 * t620 - t590 * t619) * t591;
	t605 = t596 * t588 + t593 * t590;
	t573 = t605 * t591;
	t575 = t593 * t588 - t596 * t590;
	t597 = cos(qJ(1));
	t594 = sin(qJ(1));
	t621 = qJD(1) * t594;
	t574 = -t588 * t619 - t590 * t620;
	t623 = t594 * t574;
	t546 = -t573 * t621 + t623 + (-qJD(1) * t575 - t567) * t597;
	t587 = pkin(12) + qJ(5);
	t585 = sin(t587);
	t586 = cos(t587);
	t607 = t597 * t573 - t594 * t575;
	t589 = sin(pkin(6));
	t625 = t589 * t597;
	t604 = t585 * t607 + t586 * t625;
	t612 = t589 * t621;
	t540 = t604 * qJD(5) - t546 * t586 - t585 * t612;
	t568 = qJD(2) * t573;
	t572 = t575 * t591;
	t603 = t575 * qJD(2);
	t545 = t572 * t621 + (-qJD(1) * t605 - t568) * t597 + t594 * t603;
	t613 = t585 * t625;
	t551 = -t586 * t607 + t613;
	t555 = -t597 * t572 - t594 * t605;
	t592 = sin(qJ(6));
	t595 = cos(qJ(6));
	t634 = t540 * t595 + t545 * t592 + (-t551 * t592 + t555 * t595) * qJD(6);
	t633 = (t551 * t595 + t555 * t592) * qJD(6) + t540 * t592 - t545 * t595;
	t626 = t589 * t594;
	t618 = qJD(5) * t585;
	t617 = qJD(5) * t586;
	t616 = qJD(6) * t586;
	t615 = qJD(6) * t592;
	t614 = qJD(6) * t595;
	t611 = qJD(1) * t625;
	t558 = t594 * t572 - t597 * t605;
	t598 = -t607 * qJD(1) + t594 * t567 + t597 * t574;
	t610 = -t558 * t616 + t598;
	t606 = -t594 * t573 - t597 * t575;
	t609 = t606 * qJD(1) - t555 * t616 - t597 * t567 + t623;
	t566 = t589 * t603;
	t570 = t575 * t589;
	t608 = t570 * t616 - t566;
	t571 = t605 * t589;
	t561 = t571 * t586 + t591 * t585;
	t560 = -t571 * t585 + t591 * t586;
	t552 = -t585 * t606 + t586 * t626;
	t553 = t585 * t626 + t586 * t606;
	t538 = qJD(5) * t613 - t546 * t585 + t586 * t612 - t607 * t617;
	t542 = t555 * qJD(1) - t594 * t568 - t597 * t603;
	t601 = -qJD(6) * t606 + t542 * t586 + t558 * t618;
	t600 = -qJD(6) * t607 - t545 * t586 + t555 * t618;
	t565 = qJD(2) * t571;
	t599 = qJD(6) * t571 - t565 * t586 + t570 * t618;
	t548 = t560 * qJD(5) - t566 * t586;
	t547 = -t561 * qJD(5) + t566 * t585;
	t537 = t552 * qJD(5) + t585 * t611 + t586 * t598;
	t536 = t553 * qJD(5) + t585 * t598 - t586 * t611;
	t535 = t537 * t595 + t542 * t592 + (-t553 * t592 - t558 * t595) * qJD(6);
	t534 = -t537 * t592 + t542 * t595 + (-t553 * t595 + t558 * t592) * qJD(6);
	t1 = [t634, t610 * t592 - t601 * t595, 0, 0, -t536 * t595 - t552 * t615, t534; t535, t609 * t592 - t600 * t595, 0, 0, t538 * t595 + t604 * t615, t633; 0, t608 * t592 + t599 * t595, 0, 0, t547 * t595 - t560 * t615, -t548 * t592 + t565 * t595 + (-t561 * t595 - t570 * t592) * qJD(6); -t633, t601 * t592 + t610 * t595, 0, 0, t536 * t592 - t552 * t614, -t535; t534, t600 * t592 + t609 * t595, 0, 0, -t538 * t592 + t604 * t614, t634; 0, -t599 * t592 + t608 * t595, 0, 0, -t547 * t592 - t560 * t614, -t548 * t595 - t565 * t592 + (t561 * t592 - t570 * t595) * qJD(6); t538, -t542 * t585 + t558 * t617, 0, 0, t537, 0; t536, t545 * t585 + t555 * t617, 0, 0, -t540, 0; 0, -t565 * t585 - t570 * t617, 0, 0, t548, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end