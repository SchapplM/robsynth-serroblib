% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR11
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:10
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRRPRR11_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR11_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:18
	% EndTime: 2019-10-10 12:10:18
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
	% StartTime: 2019-10-10 12:10:19
	% EndTime: 2019-10-10 12:10:19
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
	% StartTime: 2019-10-10 12:10:20
	% EndTime: 2019-10-10 12:10:20
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (94->35), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->37)
	t284 = cos(pkin(6));
	t287 = sin(qJ(1));
	t286 = sin(qJ(2));
	t307 = t287 * t286;
	t298 = t284 * t307;
	t302 = qJD(2) * t286;
	t289 = cos(qJ(2));
	t290 = cos(qJ(1));
	t304 = t290 * t289;
	t275 = -qJD(1) * t298 - t287 * t302 + (qJD(2) * t284 + qJD(1)) * t304;
	t305 = t290 * t286;
	t306 = t287 * t289;
	t277 = t284 * t305 + t306;
	t285 = sin(qJ(3));
	t288 = cos(qJ(3));
	t283 = sin(pkin(6));
	t303 = qJD(1) * t283;
	t297 = t287 * t303;
	t308 = t283 * t290;
	t311 = (-t277 * t288 + t285 * t308) * qJD(3) - t275 * t285 + t288 * t297;
	t310 = t283 * t285;
	t309 = t283 * t288;
	t301 = qJD(3) * t285;
	t300 = qJD(3) * t288;
	t299 = qJD(3) * t289;
	t296 = t290 * t303;
	t295 = t283 * qJD(2) * t289;
	t276 = t284 * t304 - t307;
	t278 = -t284 * t306 - t305;
	t293 = t298 - t304;
	t291 = -t275 * t288 + t300 * t308 + (qJD(3) * t277 - t297) * t285;
	t274 = t278 * qJD(1) - t277 * qJD(2);
	t273 = -t277 * qJD(1) + t278 * qJD(2);
	t272 = -t276 * qJD(1) + t293 * qJD(2);
	t271 = t285 * t296 + t273 * t288 + (t285 * t293 + t287 * t309) * qJD(3);
	t270 = t288 * t296 - t273 * t285 + (-t287 * t310 + t288 * t293) * qJD(3);
	t1 = [t291, t272 * t288 - t278 * t301, t270, 0, 0, 0; t271, t274 * t288 - t276 * t301, t311, 0, 0, 0; 0, (-t285 * t299 - t288 * t302) * t283, -t285 * t295 + (-t284 * t285 - t286 * t309) * qJD(3), 0, 0, 0; -t311, -t272 * t285 - t278 * t300, -t271, 0, 0, 0; t270, -t274 * t285 - t276 * t300, t291, 0, 0, 0; 0, (t285 * t302 - t288 * t299) * t283, -t288 * t295 + (-t284 * t288 + t286 * t310) * qJD(3), 0, 0, 0; t274, t273, 0, 0, 0, 0; -t272, t275, 0, 0, 0, 0; 0, t295, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:21
	% EndTime: 2019-10-10 12:10:21
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (94->35), mult. (310->70), div. (0->0), fcn. (322->8), ass. (0->37)
	t350 = cos(pkin(6));
	t353 = sin(qJ(1));
	t352 = sin(qJ(2));
	t373 = t353 * t352;
	t364 = t350 * t373;
	t368 = qJD(2) * t352;
	t355 = cos(qJ(2));
	t356 = cos(qJ(1));
	t370 = t356 * t355;
	t340 = -qJD(1) * t364 - t353 * t368 + (qJD(2) * t350 + qJD(1)) * t370;
	t371 = t356 * t352;
	t372 = t353 * t355;
	t342 = t350 * t371 + t372;
	t351 = sin(qJ(3));
	t354 = cos(qJ(3));
	t349 = sin(pkin(6));
	t369 = qJD(1) * t349;
	t363 = t353 * t369;
	t374 = t349 * t356;
	t377 = (t342 * t351 + t354 * t374) * qJD(3) - t340 * t354 - t351 * t363;
	t376 = t349 * t351;
	t375 = t349 * t354;
	t367 = qJD(3) * t351;
	t366 = qJD(3) * t354;
	t365 = qJD(3) * t355;
	t362 = t356 * t369;
	t361 = t349 * qJD(2) * t355;
	t341 = t350 * t370 - t373;
	t343 = -t350 * t372 - t371;
	t359 = t364 - t370;
	t357 = -t340 * t351 - t342 * t366 + t354 * t363 + t367 * t374;
	t339 = t343 * qJD(1) - t342 * qJD(2);
	t338 = -t342 * qJD(1) + t343 * qJD(2);
	t337 = -t341 * qJD(1) + t359 * qJD(2);
	t336 = t351 * t362 + t338 * t354 + (t351 * t359 + t353 * t375) * qJD(3);
	t335 = -t354 * t362 + t338 * t351 + (t353 * t376 - t354 * t359) * qJD(3);
	t1 = [t377, t337 * t354 - t343 * t367, -t335, 0, 0, 0; t336, t339 * t354 - t341 * t367, t357, 0, 0, 0; 0, (-t351 * t365 - t354 * t368) * t349, -t351 * t361 + (-t350 * t351 - t352 * t375) * qJD(3), 0, 0, 0; t339, t338, 0, 0, 0, 0; -t337, t340, 0, 0, 0, 0; 0, t361, 0, 0, 0, 0; t357, t337 * t351 + t343 * t366, t336, 0, 0, 0; t335, t339 * t351 + t341 * t366, -t377, 0, 0, 0; 0, (-t351 * t368 + t354 * t365) * t349, t354 * t361 + (t350 * t354 - t352 * t376) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:21
	% EndTime: 2019-10-10 12:10:21
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (371->63), mult. (1138->104), div. (0->0), fcn. (1234->10), ass. (0->55)
	t422 = cos(pkin(6));
	t426 = sin(qJ(1));
	t425 = sin(qJ(2));
	t459 = t426 * t425;
	t453 = t422 * t459;
	t454 = qJD(2) * t425;
	t429 = cos(qJ(2));
	t430 = cos(qJ(1));
	t456 = t430 * t429;
	t399 = -qJD(1) * t453 - t426 * t454 + (qJD(2) * t422 + qJD(1)) * t456;
	t457 = t430 * t425;
	t458 = t426 * t429;
	t411 = t422 * t457 + t458;
	t424 = sin(qJ(3));
	t428 = cos(qJ(3));
	t421 = sin(pkin(6));
	t460 = t421 * t430;
	t405 = -t411 * t428 + t424 * t460;
	t455 = qJD(1) * t421;
	t451 = t426 * t455;
	t392 = -t405 * qJD(3) + t399 * t424 - t428 * t451;
	t402 = t411 * t424 + t428 * t460;
	t395 = t402 * qJD(3) - t399 * t428 - t424 * t451;
	t423 = sin(qJ(5));
	t427 = cos(qJ(5));
	t465 = -t392 * t427 - t395 * t423 + (t402 * t423 - t405 * t427) * qJD(5);
	t464 = (t402 * t427 + t405 * t423) * qJD(5) - t395 * t427 + t392 * t423;
	t469 = qJD(3) - qJD(5);
	t461 = t421 * t428;
	t409 = t422 * t424 + t425 * t461;
	t449 = t421 * qJD(2) * t429;
	t400 = t409 * qJD(3) + t424 * t449;
	t462 = t421 * t424;
	t408 = -t422 * t428 + t425 * t462;
	t401 = -t408 * qJD(3) + t428 * t449;
	t468 = -t400 * t427 + t401 * t423 + (t408 * t423 + t409 * t427) * qJD(5);
	t467 = t400 * t423 + t401 * t427 + (t408 * t427 - t409 * t423) * qJD(5);
	t435 = t422 * t458 + t457;
	t397 = t411 * qJD(1) + t435 * qJD(2);
	t434 = t453 - t456;
	t407 = t426 * t462 - t428 * t434;
	t450 = t430 * t455;
	t390 = t407 * qJD(3) - t397 * t424 - t428 * t450;
	t436 = t424 * t434 + t426 * t461;
	t391 = t436 * qJD(3) - t397 * t428 + t424 * t450;
	t466 = (t407 * t427 - t423 * t436) * qJD(5) - t390 * t427 + t391 * t423;
	t438 = t423 * t428 - t424 * t427;
	t437 = t423 * t424 + t427 * t428;
	t410 = t422 * t456 - t459;
	t389 = t390 * t423 + t391 * t427 + (-t407 * t423 - t427 * t436) * qJD(5);
	t432 = t469 * t438;
	t431 = t469 * t437;
	t398 = t435 * qJD(1) + t411 * qJD(2);
	t396 = -t410 * qJD(1) + t434 * qJD(2);
	t1 = [-t464, t437 * t396 - t432 * t435, t466, 0, -t466, 0; t389, -t437 * t398 + t432 * t410, t465, 0, -t465, 0; 0, (t432 * t429 - t437 * t454) * t421, t468, 0, -t468, 0; t465, -t438 * t396 - t431 * t435, t389, 0, -t389, 0; -t466, t438 * t398 + t431 * t410, t464, 0, -t464, 0; 0, (t431 * t429 + t438 * t454) * t421, t467, 0, -t467, 0; t398, t397, 0, 0, 0, 0; t396, -t399, 0, 0, 0, 0; 0, -t449, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:10:27
	% EndTime: 2019-10-10 12:10:28
	% DurationCPUTime: 1.26s
	% Computational Cost: add. (859->107), mult. (2572->186), div. (0->0), fcn. (2868->12), ass. (0->82)
	t773 = sin(qJ(1));
	t772 = sin(qJ(2));
	t815 = cos(pkin(6));
	t795 = t773 * t815;
	t793 = t772 * t795;
	t803 = qJD(2) * t772;
	t777 = cos(qJ(2));
	t778 = cos(qJ(1));
	t805 = t778 * t777;
	t743 = -qJD(1) * t793 - t773 * t803 + (qJD(2) * t815 + qJD(1)) * t805;
	t794 = t778 * t815;
	t756 = t772 * t794 + t773 * t777;
	t771 = sin(qJ(3));
	t776 = cos(qJ(3));
	t768 = sin(pkin(6));
	t806 = t768 * t778;
	t749 = -t756 * t776 + t771 * t806;
	t804 = qJD(1) * t768;
	t799 = t773 * t804;
	t723 = -t749 * qJD(3) + t743 * t771 - t776 * t799;
	t746 = t756 * t771 + t776 * t806;
	t726 = t746 * qJD(3) - t743 * t776 - t771 * t799;
	t770 = sin(qJ(5));
	t775 = cos(qJ(5));
	t729 = t746 * t775 + t749 * t770;
	t714 = -t729 * qJD(5) - t723 * t770 + t726 * t775;
	t730 = t746 * t770 - t749 * t775;
	t782 = t778 * t772 + t777 * t795;
	t742 = t782 * qJD(1) + t756 * qJD(2);
	t755 = t773 * t772 - t777 * t794;
	t769 = sin(qJ(6));
	t774 = cos(qJ(6));
	t836 = t714 * t774 + t742 * t769 + (t730 * t769 + t755 * t774) * qJD(6);
	t835 = (t730 * t774 - t755 * t769) * qJD(6) - t714 * t769 + t742 * t774;
	t711 = t730 * qJD(5) - t723 * t775 - t726 * t770;
	t801 = qJD(6) * t774;
	t832 = t711 * t769 - t729 * t801;
	t802 = qJD(6) * t769;
	t831 = t711 * t774 + t729 * t802;
	t809 = t768 * t772;
	t753 = t771 * t809 - t815 * t776;
	t808 = t768 * t776;
	t754 = t815 * t771 + t772 * t808;
	t739 = t753 * t770 + t754 * t775;
	t807 = t768 * t777;
	t796 = qJD(2) * t807;
	t744 = t754 * qJD(3) + t771 * t796;
	t745 = -t753 * qJD(3) + t776 * t796;
	t717 = t739 * qJD(5) - t744 * t775 + t745 * t770;
	t738 = t753 * t775 - t754 * t770;
	t830 = t717 * t769 - t738 * t801;
	t829 = t717 * t774 + t738 * t802;
	t787 = t770 * t776 - t771 * t775;
	t816 = qJD(5) - qJD(3);
	t780 = t816 * t787;
	t783 = t793 - t805;
	t751 = t773 * t768 * t771 - t776 * t783;
	t785 = t771 * t783 + t773 * t808;
	t733 = -t751 * t770 - t775 * t785;
	t741 = t756 * qJD(1) + t782 * qJD(2);
	t798 = t778 * t804;
	t721 = t751 * qJD(3) - t741 * t771 - t776 * t798;
	t722 = t785 * qJD(3) - t741 * t776 + t771 * t798;
	t734 = t751 * t775 - t770 * t785;
	t781 = t734 * qJD(5) - t721 * t775 + t722 * t770;
	t820 = t733 * t801 - t769 * t781;
	t819 = t733 * t802 + t774 * t781;
	t786 = t770 * t771 + t775 * t776;
	t735 = t786 * t755;
	t710 = t733 * qJD(5) + t721 * t770 + t722 * t775;
	t719 = t738 * qJD(5) + t744 * t770 + t745 * t775;
	t797 = t768 * t803;
	t779 = t816 * t786;
	t752 = t786 * t807;
	t740 = t755 * qJD(1) + t783 * qJD(2);
	t736 = t786 * t782;
	t727 = (-t780 * t777 - t786 * t803) * t768;
	t716 = -t786 * t742 + t755 * t780;
	t715 = t786 * t740 + t780 * t782;
	t707 = t710 * t774 + t740 * t769 + (-t734 * t769 - t774 * t782) * qJD(6);
	t706 = -t710 * t769 + t740 * t774 + (-t734 * t774 + t769 * t782) * qJD(6);
	t1 = [t836, t715 * t774 + t741 * t769 + (t736 * t769 + t774 * t783) * qJD(6), t819, 0, -t819, t706; t707, t716 * t774 - t743 * t769 + (t735 * t769 - t756 * t774) * qJD(6), t831, 0, -t831, -t835; 0, -t769 * t796 + t727 * t774 + (-t752 * t769 - t774 * t809) * qJD(6), t829, 0, -t829, -t774 * t797 - t719 * t769 + (-t739 * t774 - t769 * t807) * qJD(6); t835, -t715 * t769 + t741 * t774 + (t736 * t774 - t769 * t783) * qJD(6), t820, 0, -t820, -t707; t706, -t716 * t769 - t743 * t774 + (t735 * t774 + t756 * t769) * qJD(6), -t832, 0, t832, t836; 0, -t774 * t796 - t727 * t769 + (-t752 * t774 + t769 * t809) * qJD(6), -t830, 0, t830, t769 * t797 - t719 * t774 + (t739 * t769 - t774 * t807) * qJD(6); -t711, t787 * t740 - t779 * t782, -t710, 0, t710, 0; t781, -t816 * t735 - t787 * t742, t714, 0, -t714, 0; 0, (t779 * t777 - t787 * t803) * t768, -t719, 0, t719, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end