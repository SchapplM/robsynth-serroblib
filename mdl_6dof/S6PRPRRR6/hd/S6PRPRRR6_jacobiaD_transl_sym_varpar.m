% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRRR6
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:03
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(11));
	t50 = sin(pkin(11));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (24->13), mult. (82->26), div. (0->0), fcn. (72->6), ass. (0->16)
	t135 = -pkin(2) + r_i_i_C(2);
	t134 = r_i_i_C(3) + qJ(3);
	t126 = cos(pkin(6));
	t127 = sin(qJ(2));
	t133 = t126 * t127;
	t128 = cos(qJ(2));
	t132 = t126 * t128;
	t131 = qJD(2) * t134;
	t123 = sin(pkin(11));
	t125 = cos(pkin(11));
	t130 = t123 * t128 + t125 * t133;
	t129 = t123 * t133 - t125 * t128;
	t124 = sin(pkin(6));
	t122 = t129 * qJD(2);
	t120 = t130 * qJD(2);
	t1 = [0, -t129 * qJD(3) - t135 * t122 - (t123 * t132 + t125 * t127) * t131, -t122, 0, 0, 0; 0, t130 * qJD(3) + t135 * t120 - (t123 * t127 - t125 * t132) * t131, t120, 0, 0, 0; 0, (qJD(3) * t127 + (t135 * t127 + t134 * t128) * qJD(2)) * t124, t124 * qJD(2) * t127, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:07
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (67->29), mult. (228->61), div. (0->0), fcn. (208->8), ass. (0->28)
	t178 = sin(pkin(6));
	t181 = sin(qJ(4));
	t199 = t178 * t181;
	t183 = cos(qJ(4));
	t198 = t178 * t183;
	t184 = cos(qJ(2));
	t197 = t178 * t184;
	t179 = cos(pkin(11));
	t196 = t179 * t184;
	t180 = cos(pkin(6));
	t182 = sin(qJ(2));
	t195 = t180 * t182;
	t194 = t180 * t184;
	t193 = qJD(2) * t182;
	t192 = -pkin(2) - pkin(8) - r_i_i_C(3);
	t177 = sin(pkin(11));
	t191 = t177 * t193;
	t190 = t178 * t193;
	t189 = qJD(2) * t196;
	t188 = t183 * r_i_i_C(1) - t181 * r_i_i_C(2);
	t187 = t181 * r_i_i_C(1) + t183 * r_i_i_C(2) + qJ(3);
	t186 = t177 * t184 + t179 * t195;
	t173 = t177 * t194 + t179 * t182;
	t185 = t188 * qJD(4) + qJD(3);
	t171 = t177 * t182 - t179 * t194;
	t170 = -t180 * t191 + t189;
	t168 = t186 * qJD(2);
	t1 = [0, t185 * (-t177 * t195 + t196) + t192 * t170 - t187 * t173 * qJD(2), t170, t188 * t170 + ((-t173 * t181 - t177 * t198) * r_i_i_C(1) + (-t173 * t183 + t177 * t199) * r_i_i_C(2)) * qJD(4), 0, 0; 0, t185 * t186 + t192 * t168 - t187 * (-t180 * t189 + t191), t168, t188 * t168 + ((-t171 * t181 + t179 * t198) * r_i_i_C(1) + (-t171 * t183 - t179 * t199) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (t185 * t182 + (t192 * t182 + t187 * t184) * qJD(2)) * t178, t190, t188 * t190 + ((-t180 * t183 + t181 * t197) * r_i_i_C(1) + (t180 * t181 + t183 * t197) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:08
	% EndTime: 2019-10-09 22:03:09
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (236->64), mult. (762->121), div. (0->0), fcn. (760->10), ass. (0->49)
	t293 = sin(pkin(11));
	t295 = cos(pkin(11));
	t299 = sin(qJ(2));
	t296 = cos(pkin(6));
	t302 = cos(qJ(2));
	t322 = t296 * t302;
	t285 = t293 * t322 + t295 * t299;
	t301 = cos(qJ(4));
	t294 = sin(pkin(6));
	t298 = sin(qJ(4));
	t327 = t294 * t298;
	t331 = -t285 * t301 + t293 * t327;
	t297 = sin(qJ(5));
	t300 = cos(qJ(5));
	t314 = r_i_i_C(1) * t300 - r_i_i_C(2) * t297;
	t306 = t314 * qJD(5);
	t312 = pkin(4) + t314;
	t329 = pkin(9) + r_i_i_C(3);
	t330 = t312 * t298 - t329 * t301 + qJ(3);
	t326 = t294 * t299;
	t325 = t294 * t301;
	t324 = t294 * t302;
	t323 = t296 * t299;
	t321 = qJD(2) * t299;
	t320 = qJD(2) * t302;
	t318 = t293 * t321;
	t317 = t294 * t320;
	t316 = t295 * t320;
	t315 = t294 * t321;
	t313 = -r_i_i_C(1) * t297 - r_i_i_C(2) * t300;
	t311 = -pkin(2) - pkin(8) + t313;
	t283 = t293 * t299 - t295 * t322;
	t310 = -t283 * t298 + t295 * t325;
	t309 = t283 * t301 + t295 * t327;
	t272 = t285 * t298 + t293 * t325;
	t284 = t293 * t302 + t295 * t323;
	t308 = t296 * t298 + t301 * t324;
	t307 = -t296 * t301 + t298 * t324;
	t305 = qJD(5) * t313;
	t303 = qJD(3) + t298 * t305 + (t329 * t298 + t312 * t301) * qJD(4);
	t286 = -t293 * t323 + t295 * t302;
	t282 = -t296 * t318 + t316;
	t281 = t285 * qJD(2);
	t280 = t284 * qJD(2);
	t279 = -t296 * t316 + t318;
	t275 = t308 * qJD(4) - t298 * t315;
	t269 = t309 * qJD(4) + t280 * t298;
	t267 = t331 * qJD(4) - t282 * t298;
	t1 = [0, -t281 * t330 + t311 * t282 - t285 * t306 + t303 * t286, t282, -t329 * t267 - t331 * t305 + t312 * (-t272 * qJD(4) + t282 * t301), (t267 * t297 - t281 * t300) * r_i_i_C(1) + (t267 * t300 + t281 * t297) * r_i_i_C(2) + ((-t272 * t300 - t286 * t297) * r_i_i_C(1) + (t272 * t297 - t286 * t300) * r_i_i_C(2)) * qJD(5), 0; 0, -t279 * t330 + t311 * t280 - t283 * t306 + t303 * t284, t280, t329 * t269 + t309 * t305 + t312 * (t310 * qJD(4) + t280 * t301), (-t269 * t297 - t279 * t300) * r_i_i_C(1) + (-t269 * t300 + t279 * t297) * r_i_i_C(2) + ((-t284 * t297 + t300 * t310) * r_i_i_C(1) + (-t284 * t300 - t297 * t310) * r_i_i_C(2)) * qJD(5), 0; 0, ((t330 * qJD(2) + t306) * t302 + (t311 * qJD(2) + t303) * t299) * t294, t315, -t329 * t275 - t308 * t305 + t312 * (t307 * qJD(4) + t301 * t315), (t275 * t297 + t300 * t317) * r_i_i_C(1) + (t275 * t300 - t297 * t317) * r_i_i_C(2) + ((-t297 * t326 + t300 * t307) * r_i_i_C(1) + (-t297 * t307 - t300 * t326) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:08
	% EndTime: 2019-10-09 22:03:09
	% DurationCPUTime: 0.55s
	% Computational Cost: add. (499->77), mult. (1130->129), div. (0->0), fcn. (1134->12), ass. (0->63)
	t335 = sin(pkin(11));
	t337 = cos(pkin(11));
	t340 = sin(qJ(2));
	t343 = cos(qJ(2));
	t381 = cos(pkin(6));
	t365 = t343 * t381;
	t322 = t335 * t365 + t337 * t340;
	t342 = cos(qJ(4));
	t336 = sin(pkin(6));
	t339 = sin(qJ(4));
	t379 = t336 * t339;
	t385 = -t322 * t342 + t335 * t379;
	t333 = qJD(5) + qJD(6);
	t341 = cos(qJ(5));
	t334 = qJ(5) + qJ(6);
	t331 = sin(t334);
	t332 = cos(t334);
	t358 = r_i_i_C(1) * t332 - r_i_i_C(2) * t331;
	t382 = pkin(5) * qJD(5);
	t347 = -t358 * t333 - t341 * t382;
	t355 = pkin(5) * t341 + pkin(4) + t358;
	t383 = r_i_i_C(3) + pkin(10) + pkin(9);
	t384 = t339 * t355 - t342 * t383 + qJ(3);
	t378 = t336 * t340;
	t377 = t336 * t342;
	t376 = t336 * t343;
	t309 = t322 * t339 + t335 * t377;
	t318 = t322 * qJD(2);
	t362 = t309 * t333 + t318;
	t366 = t340 * t381;
	t360 = t335 * t366;
	t371 = qJD(2) * t343;
	t319 = -qJD(2) * t360 + t337 * t371;
	t304 = qJD(4) * t385 - t319 * t339;
	t323 = t337 * t343 - t360;
	t364 = -t323 * t333 + t304;
	t375 = (t331 * t364 - t332 * t362) * r_i_i_C(1) + (t331 * t362 + t332 * t364) * r_i_i_C(2);
	t359 = t337 * t365;
	t372 = qJD(2) * t340;
	t316 = -qJD(2) * t359 + t335 * t372;
	t320 = t335 * t340 - t359;
	t354 = -t320 * t339 + t337 * t377;
	t361 = -t333 * t354 + t316;
	t321 = t335 * t343 + t337 * t366;
	t317 = t321 * qJD(2);
	t353 = t320 * t342 + t337 * t379;
	t306 = qJD(4) * t353 + t317 * t339;
	t363 = -t321 * t333 - t306;
	t374 = (t331 * t363 - t332 * t361) * r_i_i_C(1) + (t331 * t361 + t332 * t363) * r_i_i_C(2);
	t351 = t339 * t376 - t342 * t381;
	t368 = t336 * t371;
	t352 = t333 * t351 + t368;
	t350 = t339 * t381 + t342 * t376;
	t367 = t336 * t372;
	t312 = qJD(4) * t350 - t339 * t367;
	t356 = -t333 * t378 + t312;
	t373 = (t331 * t356 + t332 * t352) * r_i_i_C(1) + (-t331 * t352 + t332 * t356) * r_i_i_C(2);
	t357 = -r_i_i_C(1) * t331 - r_i_i_C(2) * t332;
	t338 = sin(qJ(5));
	t349 = -pkin(5) * t338 - pkin(2) - pkin(8) + t357;
	t348 = t333 * t357 - t338 * t382;
	t345 = qJD(3) + t348 * t339 + (t339 * t383 + t342 * t355) * qJD(4);
	t1 = [0, -t318 * t384 + t319 * t349 + t322 * t347 + t323 * t345, t319, -t383 * t304 - t348 * t385 + t355 * (-qJD(4) * t309 + t319 * t342), (t304 * t338 - t318 * t341 + (-t309 * t341 - t323 * t338) * qJD(5)) * pkin(5) + t375, t375; 0, -t316 * t384 + t317 * t349 + t320 * t347 + t321 * t345, t317, t383 * t306 + t348 * t353 + t355 * (qJD(4) * t354 + t317 * t342), (-t306 * t338 - t316 * t341 + (-t321 * t338 + t341 * t354) * qJD(5)) * pkin(5) + t374, t374; 0, ((qJD(2) * t384 - t347) * t343 + (qJD(2) * t349 + t345) * t340) * t336, t367, -t383 * t312 - t348 * t350 + t355 * (qJD(4) * t351 + t342 * t367), (t341 * t368 + t312 * t338 + (-t338 * t378 + t341 * t351) * qJD(5)) * pkin(5) + t373, t373;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end