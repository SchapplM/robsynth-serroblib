% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRPR10_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR10_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
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
	% StartTime: 2019-10-10 10:20:53
	% EndTime: 2019-10-10 10:20:53
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
	% StartTime: 2019-10-10 10:20:53
	% EndTime: 2019-10-10 10:20:53
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
	t229 = cos(pkin(11));
	t227 = sin(pkin(11));
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
	% StartTime: 2019-10-10 10:20:54
	% EndTime: 2019-10-10 10:20:54
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
	t301 = pkin(11) + qJ(4);
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
	% StartTime: 2019-10-10 10:20:55
	% EndTime: 2019-10-10 10:20:55
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (144->36), mult. (310->71), div. (0->0), fcn. (322->8), ass. (0->38)
	t360 = cos(pkin(6));
	t362 = sin(qJ(1));
	t361 = sin(qJ(2));
	t381 = t362 * t361;
	t372 = t360 * t381;
	t376 = qJD(2) * t361;
	t363 = cos(qJ(2));
	t364 = cos(qJ(1));
	t378 = t364 * t363;
	t348 = -qJD(1) * t372 - t362 * t376 + (qJD(2) * t360 + qJD(1)) * t378;
	t379 = t364 * t361;
	t380 = t362 * t363;
	t350 = t360 * t379 + t380;
	t358 = pkin(11) + qJ(4);
	t356 = sin(t358);
	t357 = cos(t358);
	t359 = sin(pkin(6));
	t377 = qJD(1) * t359;
	t371 = t362 * t377;
	t382 = t359 * t364;
	t385 = (-t350 * t357 + t356 * t382) * qJD(4) - t348 * t356 + t357 * t371;
	t384 = t359 * t361;
	t383 = t359 * t362;
	t375 = qJD(4) * t356;
	t374 = qJD(4) * t357;
	t373 = qJD(4) * t363;
	t370 = t364 * t377;
	t369 = t359 * qJD(2) * t363;
	t349 = t360 * t378 - t381;
	t351 = -t360 * t380 - t379;
	t367 = t372 - t378;
	t365 = t348 * t357 + t356 * t371 + (-t350 * t356 - t357 * t382) * qJD(4);
	t347 = t351 * qJD(1) - t350 * qJD(2);
	t346 = -t350 * qJD(1) + t351 * qJD(2);
	t345 = -t349 * qJD(1) + t367 * qJD(2);
	t344 = t356 * t370 + t346 * t357 + (t356 * t367 + t357 * t383) * qJD(4);
	t343 = -t357 * t370 + t346 * t356 + (t356 * t383 - t357 * t367) * qJD(4);
	t1 = [t347, t346, 0, 0, 0, 0; -t345, t348, 0, 0, 0, 0; 0, t369, 0, 0, 0, 0; t365, -t345 * t357 + t351 * t375, 0, t343, 0, 0; -t344, -t347 * t357 + t349 * t375, 0, -t385, 0, 0; 0, (t356 * t373 + t357 * t376) * t359, 0, t356 * t369 + (t356 * t360 + t357 * t384) * qJD(4), 0, 0; t385, t345 * t356 + t351 * t374, 0, t344, 0, 0; t343, t347 * t356 + t349 * t374, 0, t365, 0, 0; 0, (-t356 * t376 + t357 * t373) * t359, 0, t357 * t369 + (-t356 * t384 + t357 * t360) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:56
	% EndTime: 2019-10-10 10:20:56
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (424->76), mult. (898->140), div. (0->0), fcn. (964->10), ass. (0->68)
	t479 = sin(qJ(1));
	t476 = cos(pkin(6));
	t491 = qJD(2) * t476 + qJD(1);
	t478 = sin(qJ(2));
	t512 = t479 * t478;
	t498 = t476 * t512;
	t506 = qJD(2) * t478;
	t481 = cos(qJ(2));
	t482 = cos(qJ(1));
	t508 = t482 * t481;
	t451 = -qJD(1) * t498 - t479 * t506 + t491 * t508;
	t509 = t482 * t478;
	t511 = t479 * t481;
	t462 = t476 * t509 + t511;
	t474 = pkin(11) + qJ(4);
	t472 = sin(t474);
	t473 = cos(t474);
	t475 = sin(pkin(6));
	t507 = qJD(1) * t475;
	t495 = t479 * t507;
	t514 = t475 * t482;
	t497 = t472 * t514;
	t503 = qJD(4) * t473;
	t445 = -qJD(4) * t497 + t451 * t472 + t462 * t503 - t473 * t495;
	t463 = t476 * t511 + t509;
	t450 = t463 * qJD(1) + t462 * qJD(2);
	t454 = t462 * t472 + t473 * t514;
	t496 = t476 * t508;
	t461 = -t496 + t512;
	t477 = sin(qJ(6));
	t480 = cos(qJ(6));
	t525 = (t454 * t477 + t461 * t480) * qJD(6) - t445 * t480 + t450 * t477;
	t522 = -t445 * t477 - t450 * t480 + (-t454 * t480 + t461 * t477) * qJD(6);
	t521 = t454 * qJD(4) - t451 * t473 - t472 * t495;
	t516 = t475 * t478;
	t515 = t475 * t479;
	t513 = t477 * t481;
	t510 = t480 * t481;
	t505 = qJD(2) * t481;
	t504 = qJD(4) * t472;
	t502 = qJD(4) * t481;
	t501 = qJD(6) * t472;
	t500 = qJD(6) * t477;
	t499 = qJD(6) * t480;
	t494 = t482 * t507;
	t493 = t475 * t506;
	t492 = t475 * t505;
	t490 = qJD(2) + t501;
	t449 = -t462 * qJD(1) - t463 * qJD(2);
	t489 = t463 * t501 - t449;
	t488 = t461 * t501 - t451;
	t464 = -t498 + t508;
	t487 = -t464 * t472 + t473 * t515;
	t458 = t464 * t473 + t472 * t515;
	t460 = t476 * t472 + t473 * t516;
	t459 = t472 * t516 - t476 * t473;
	t448 = -qJD(1) * t496 - t482 * t505 + t491 * t512;
	t485 = -qJD(6) * t464 + t448 * t472 - t463 * t503;
	t484 = -qJD(6) * t462 - t450 * t472 - t461 * t503;
	t483 = t473 * t502 + (-qJD(2) * t472 - qJD(6)) * t478;
	t455 = t462 * t473 - t497;
	t453 = -t459 * qJD(4) + t473 * t492;
	t452 = t460 * qJD(4) + t472 * t492;
	t444 = t487 * qJD(4) + t449 * t473 + t472 * t494;
	t443 = t458 * qJD(4) + t449 * t472 - t473 * t494;
	t442 = t443 * t477 - t448 * t480 + (-t463 * t477 - t480 * t487) * qJD(6);
	t441 = t443 * t480 + t448 * t477 + (-t463 * t480 + t477 * t487) * qJD(6);
	t1 = [t522, t485 * t477 - t489 * t480, 0, t444 * t477 + t458 * t499, 0, t441; t442, t484 * t477 - t488 * t480, 0, t455 * t499 - t477 * t521, 0, -t525; 0, (t483 * t477 + t490 * t510) * t475, 0, t453 * t477 + t460 * t499, 0, -t477 * t493 + t452 * t480 + (-t459 * t477 + t475 * t510) * qJD(6); t525, t489 * t477 + t485 * t480, 0, t444 * t480 - t458 * t500, 0, -t442; t441, t488 * t477 + t484 * t480, 0, -t455 * t500 - t480 * t521, 0, t522; 0, (t483 * t480 - t490 * t513) * t475, 0, t453 * t480 - t460 * t500, 0, -t480 * t493 - t452 * t477 + (-t459 * t480 - t475 * t513) * qJD(6); t521, t448 * t473 + t463 * t504, 0, -t443, 0, 0; t444, -t450 * t473 + t461 * t504, 0, -t445, 0, 0; 0, (-t472 * t502 - t473 * t506) * t475, 0, -t452, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end