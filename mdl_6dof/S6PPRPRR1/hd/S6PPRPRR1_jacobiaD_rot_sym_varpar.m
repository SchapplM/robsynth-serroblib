% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6PPRPRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:08
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRPRR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRPRR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobiaD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (151->10), mult. (481->28), div. (18->4), fcn. (595->10), ass. (0->20)
	t84 = sin(pkin(12));
	t88 = cos(pkin(12));
	t89 = cos(pkin(11));
	t85 = sin(pkin(11));
	t98 = t85 * cos(pkin(6));
	t82 = -t84 * t98 + t89 * t88;
	t92 = sin(qJ(3));
	t93 = cos(qJ(3));
	t96 = t85 * sin(pkin(7)) * sin(pkin(6)) + (-t89 * t84 - t88 * t98) * cos(pkin(7));
	t78 = t82 * t92 - t96 * t93;
	t79 = t82 * t93 + t96 * t92;
	t76 = 0.1e1 / t79 ^ 2;
	t104 = qJD(3) * t76;
	t101 = t79 * t104;
	t102 = t78 / t79 * t104;
	t75 = t78 ^ 2;
	t72 = t75 * t76 + 0.1e1;
	t103 = (t78 * t101 + t75 * t102) / t72 ^ 2;
	t70 = 0.1e1 / t72;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, -0.2e1 * t103 + 0.2e1 * (t70 * t101 + (t70 * t102 - t76 * t103) * t78) * t78, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:44
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (337->14), mult. (1021->37), div. (25->4), fcn. (1288->12), ass. (0->26)
	t112 = sin(pkin(12));
	t117 = cos(pkin(12));
	t118 = cos(pkin(11));
	t113 = sin(pkin(11));
	t129 = t113 * cos(pkin(6));
	t108 = -t112 * t129 + t118 * t117;
	t111 = sin(pkin(13));
	t116 = cos(pkin(13));
	t121 = sin(qJ(3));
	t122 = cos(qJ(3));
	t110 = t121 * t111 - t122 * t116;
	t127 = t122 * t111 + t121 * t116;
	t138 = sin(pkin(7)) * t113 * sin(pkin(6)) + (-t118 * t112 - t117 * t129) * cos(pkin(7));
	t125 = -t108 * t110 + t138 * t127;
	t97 = 0.1e1 / t125 ^ 2;
	t139 = t125 * t97;
	t100 = -t108 * t127 - t110 * t138;
	t96 = 0.1e1 / t125;
	t94 = t100 * qJD(3);
	t133 = t94 * t96 * t97;
	t134 = qJD(3) * t139;
	t95 = t100 ^ 2;
	t92 = t95 * t97 + 0.1e1;
	t135 = (-t100 * t134 - t95 * t133) / t92 ^ 2;
	t90 = 0.1e1 / t92;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, -0.2e1 * t125 * t96 * t135 + 0.2e1 * (-t90 * t134 + (-t90 * t133 - t97 * t135) * t100) * t100 + (t96 - t139) * t90 * t94, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:44
	% EndTime: 2019-10-09 21:08:45
	% DurationCPUTime: 1.24s
	% Computational Cost: add. (4637->79), mult. (13855->173), div. (281->12), fcn. (18065->17), ass. (0->88)
	t255 = sin(pkin(13));
	t260 = cos(pkin(13));
	t266 = sin(qJ(3));
	t268 = cos(qJ(3));
	t281 = t268 * t255 + t266 * t260;
	t298 = qJD(3) * t281;
	t251 = t266 * t255 - t268 * t260;
	t263 = cos(pkin(7));
	t243 = t251 * t263;
	t256 = sin(pkin(12));
	t257 = sin(pkin(11));
	t261 = cos(pkin(12));
	t262 = cos(pkin(11));
	t264 = cos(pkin(6));
	t286 = t262 * t264;
	t245 = -t257 * t256 + t261 * t286;
	t246 = t256 * t286 + t257 * t261;
	t259 = sin(pkin(6));
	t258 = sin(pkin(7));
	t278 = t251 * t258;
	t276 = t259 * t278;
	t221 = -t245 * t243 - t246 * t281 + t262 * t276;
	t232 = (t261 * t243 + t256 * t281) * t259 + t264 * t278;
	t207 = atan2(t221, t232);
	t202 = sin(t207);
	t203 = cos(t207);
	t196 = t202 * t221 + t203 * t232;
	t193 = 0.1e1 / t196;
	t288 = t257 * t264;
	t247 = -t262 * t256 - t261 * t288;
	t289 = t257 * t259;
	t234 = -t247 * t258 + t263 * t289;
	t265 = sin(qJ(5));
	t267 = cos(qJ(5));
	t242 = t281 * t258;
	t244 = t281 * t263;
	t248 = -t256 * t288 + t262 * t261;
	t275 = t242 * t289 + t247 * t244 - t248 * t251;
	t213 = t234 * t265 + t267 * t275;
	t209 = 0.1e1 / t213;
	t228 = 0.1e1 / t232;
	t194 = 0.1e1 / t196 ^ 2;
	t210 = 0.1e1 / t213 ^ 2;
	t229 = 0.1e1 / t232 ^ 2;
	t218 = t221 ^ 2;
	t206 = t218 * t229 + 0.1e1;
	t204 = 0.1e1 / t206;
	t239 = t258 * t298;
	t241 = t263 * t298;
	t249 = t251 * qJD(3);
	t287 = t259 * t262;
	t214 = t239 * t287 - t245 * t241 + t246 * t249;
	t226 = t264 * t239 + (t241 * t261 - t249 * t256) * t259;
	t292 = t221 * t229;
	t187 = (t214 * t228 - t226 * t292) * t204;
	t282 = -t202 * t232 + t203 * t221;
	t184 = t282 * t187 + t202 * t214 + t203 * t226;
	t297 = t184 * t193 * t194;
	t212 = -t234 * t267 + t265 * t275;
	t208 = t212 ^ 2;
	t199 = t208 * t210 + 0.1e1;
	t238 = t258 * t249;
	t240 = t263 * t249;
	t217 = -t238 * t289 - t247 * t240 - t248 * t298;
	t200 = t213 * qJD(5) + t217 * t265;
	t293 = t210 * t212;
	t283 = qJD(5) * t212;
	t201 = t217 * t267 - t283;
	t294 = t201 * t209 * t210;
	t296 = (t200 * t293 - t208 * t294) / t199 ^ 2;
	t223 = -t247 * t243 - t248 * t281 - t257 * t276;
	t295 = t194 * t223;
	t231 = t264 * t242 + (t244 * t261 - t251 * t256) * t259;
	t291 = t221 * t231;
	t290 = t226 * t228 * t229;
	t280 = -t209 * t265 + t267 * t293;
	t220 = t242 * t287 - t245 * t244 + t246 * t251;
	t279 = -t220 * t228 + t229 * t291;
	t227 = -t264 * t238 + (-t240 * t261 - t256 * t298) * t259;
	t219 = t223 ^ 2;
	t216 = -t239 * t289 - t247 * t241 + t248 * t249;
	t215 = -t238 * t287 + t245 * t240 + t246 * t298;
	t197 = 0.1e1 / t199;
	t191 = t219 * t194 + 0.1e1;
	t188 = t279 * t204;
	t185 = -t282 * t188 + t202 * t220 + t203 * t231;
	t183 = 0.2e1 * t279 / t206 ^ 2 * (t214 * t292 - t218 * t290) + (0.2e1 * t290 * t291 + t215 * t228 + (-t214 * t231 - t220 * t226 - t221 * t227) * t229) * t204;
	t1 = [0, 0, t183, 0, 0, 0; 0, 0, 0.2e1 * (-t185 * t295 - t193 * t275) / t191 ^ 2 * (t216 * t295 - t219 * t297) + (t217 * t193 + (-t184 * t275 + t185 * t216) * t194 + (-0.2e1 * t185 * t297 + ((t183 * t221 - t188 * t214 + t227 + (t188 * t232 + t220) * t187) * t203 + (-t183 * t232 + t188 * t226 + t215 + (t188 * t221 - t231) * t187) * t202) * t194) * t223) / t191, 0, 0, 0; 0, 0, 0.2e1 * t280 * t223 * t296 + (-t280 * t216 + ((qJD(5) * t209 + 0.2e1 * t212 * t294) * t267 + (-t200 * t267 + (-t201 + t283) * t265) * t210) * t223) * t197, 0, -0.2e1 * t296 + 0.2e1 * (t197 * t200 * t210 + (-t197 * t294 - t210 * t296) * t212) * t212, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:44
	% EndTime: 2019-10-09 21:08:47
	% DurationCPUTime: 3.24s
	% Computational Cost: add. (14279->140), mult. (41349->282), div. (559->12), fcn. (54415->19), ass. (0->130)
	t357 = sin(pkin(11));
	t360 = cos(pkin(12));
	t361 = cos(pkin(6));
	t356 = sin(pkin(12));
	t420 = cos(pkin(11));
	t393 = t420 * t356;
	t345 = t357 * t360 + t361 * t393;
	t355 = sin(pkin(13));
	t359 = cos(pkin(13));
	t364 = sin(qJ(3));
	t367 = cos(qJ(3));
	t350 = t355 * t364 - t367 * t359;
	t419 = sin(pkin(7));
	t394 = t367 * t419;
	t396 = t364 * t419;
	t375 = -t355 * t394 - t359 * t396;
	t421 = cos(pkin(7));
	t395 = t367 * t421;
	t397 = t364 * t421;
	t376 = -t355 * t395 - t359 * t397;
	t392 = t420 * t360;
	t381 = t357 * t356 - t361 * t392;
	t358 = sin(pkin(6));
	t398 = t358 * t420;
	t321 = -t345 * t350 + t375 * t398 + t376 * t381;
	t363 = sin(qJ(5));
	t366 = cos(qJ(5));
	t374 = t381 * t419 - t421 * t398;
	t311 = t321 * t366 + t374 * t363;
	t340 = -t355 * t396 + t359 * t394;
	t336 = t340 * qJD(3);
	t342 = -t355 * t397 + t359 * t395;
	t338 = t342 * qJD(3);
	t351 = -t355 * t367 - t364 * t359;
	t349 = t351 * qJD(3);
	t315 = -t336 * t398 - t381 * t338 + t345 * t349;
	t287 = t311 * qJD(5) + t315 * t363;
	t309 = t321 * t363 - t374 * t366;
	t307 = t309 ^ 2;
	t333 = -t361 * t375 + (-t350 * t356 - t360 * t376) * t358;
	t344 = -t358 * t360 * t419 + t361 * t421;
	t328 = t333 * t363 - t344 * t366;
	t326 = 0.1e1 / t328 ^ 2;
	t301 = t307 * t326 + 0.1e1;
	t299 = 0.1e1 / t301;
	t329 = t333 * t366 + t344 * t363;
	t331 = t361 * t336 + (t338 * t360 + t349 * t356) * t358;
	t305 = t329 * qJD(5) + t331 * t363;
	t325 = 0.1e1 / t328;
	t409 = t309 * t326;
	t271 = (-t287 * t325 + t305 * t409) * t299;
	t302 = atan2(-t309, t328);
	t297 = sin(t302);
	t298 = cos(t302);
	t385 = -t297 * t328 - t298 * t309;
	t267 = t385 * t271 - t297 * t287 + t298 * t305;
	t281 = -t297 * t309 + t298 * t328;
	t278 = 0.1e1 / t281;
	t279 = 0.1e1 / t281 ^ 2;
	t425 = t267 * t278 * t279;
	t405 = t357 * t361;
	t346 = -t360 * t405 - t393;
	t406 = t357 * t358;
	t378 = -t346 * t419 + t421 * t406;
	t347 = -t356 * t405 + t392;
	t379 = -t346 * t376 - t347 * t350 - t375 * t406;
	t312 = t363 * t379 - t378 * t366;
	t424 = 0.2e1 * t312 * t425;
	t320 = -t340 * t398 - t381 * t342 + t345 * t351;
	t332 = t361 * t340 + (t342 * t360 + t351 * t356) * t358;
	t382 = -t320 * t325 + t332 * t409;
	t423 = t363 * t382;
	t410 = t305 * t325 * t326;
	t422 = -0.2e1 * (t287 * t409 - t307 * t410) / t301 ^ 2;
	t313 = t378 * t363 + t366 * t379;
	t323 = t340 * t406 + t342 * t346 + t347 * t351;
	t362 = sin(qJ(6));
	t365 = cos(qJ(6));
	t296 = t313 * t365 - t323 * t362;
	t292 = 0.1e1 / t296;
	t293 = 0.1e1 / t296 ^ 2;
	t418 = t279 * t312;
	t380 = t336 * t406 + t338 * t346 + t347 * t349;
	t290 = -t312 * qJD(5) + t366 * t380;
	t337 = t375 * qJD(3);
	t339 = t376 * qJD(3);
	t348 = t350 * qJD(3);
	t316 = t337 * t406 + t339 * t346 + t347 * t348;
	t295 = t313 * t362 + t323 * t365;
	t403 = qJD(6) * t295;
	t283 = t290 * t365 - t316 * t362 - t403;
	t417 = t283 * t292 * t293;
	t282 = t296 * qJD(6) + t290 * t362 + t316 * t365;
	t291 = t295 ^ 2;
	t286 = t291 * t293 + 0.1e1;
	t414 = t293 * t295;
	t416 = 0.1e1 / t286 ^ 2 * (t282 * t414 - t291 * t417);
	t415 = t292 * t362;
	t413 = t295 * t365;
	t412 = t297 * t312;
	t411 = t298 * t312;
	t408 = t323 * t363;
	t407 = t323 * t366;
	t404 = qJD(5) * t366;
	t308 = t312 ^ 2;
	t277 = t279 * t308 + 0.1e1;
	t289 = t313 * qJD(5) + t363 * t380;
	t402 = 0.2e1 * (t289 * t418 - t308 * t425) / t277 ^ 2;
	t400 = -0.2e1 * t416;
	t399 = t295 * t417;
	t391 = -0.2e1 * t309 * t410;
	t386 = qJD(6) * t407 - t380;
	t384 = t293 * t413 - t415;
	t383 = -t311 * t325 + t329 * t409;
	t377 = -qJD(5) * t408 + qJD(6) * t379 + t316 * t366;
	t330 = t361 * t337 + (t339 * t360 + t348 * t356) * t358;
	t314 = -t337 * t398 - t381 * t339 + t345 * t348;
	t306 = -t328 * qJD(5) + t331 * t366;
	t304 = t362 * t379 + t365 * t407;
	t303 = t362 * t407 - t365 * t379;
	t288 = -t309 * qJD(5) + t315 * t366;
	t284 = 0.1e1 / t286;
	t275 = 0.1e1 / t277;
	t273 = t299 * t423;
	t272 = t383 * t299;
	t269 = (-t297 * t320 + t298 * t332) * t363 + t385 * t273;
	t268 = t385 * t272 - t297 * t311 + t298 * t329;
	t265 = t383 * t422 + (t329 * t391 - t288 * t325 + (t287 * t329 + t305 * t311 + t306 * t309) * t326) * t299;
	t264 = t422 * t423 + (t382 * t404 + (t332 * t391 - t314 * t325 + (t287 * t332 + t305 * t320 + t309 * t330) * t326) * t363) * t299;
	t1 = [0, 0, t264, 0, t265, 0; 0, 0, (t269 * t418 - t278 * t408) * t402 + ((t316 * t363 + t323 * t404) * t278 + t269 * t424 + (-t269 * t289 - t408 * t267 - (t332 * t404 - t264 * t309 - t273 * t287 + t330 * t363 + (-t273 * t328 - t320 * t363) * t271) * t411 - (-t320 * t404 - t264 * t328 - t273 * t305 - t314 * t363 + (t273 * t309 - t332 * t363) * t271) * t412) * t279) * t275, 0, (t268 * t418 - t278 * t313) * t402 + (t268 * t424 + t290 * t278 + (-t313 * t267 - t268 * t289 - (-t265 * t309 - t272 * t287 + t306 + (-t272 * t328 - t311) * t271) * t411 - (-t265 * t328 - t272 * t305 - t288 + (t272 * t309 - t329) * t271) * t412) * t279) * t275, 0; 0, 0, 0.2e1 * (-t292 * t303 + t304 * t414) * t416 + (0.2e1 * t304 * t399 + t386 * t292 * t365 + t377 * t415 + (t386 * t295 * t362 - t304 * t282 - t303 * t283 - t377 * t413) * t293) * t284, 0, t384 * t312 * t400 + (t384 * t289 + ((-qJD(6) * t292 - 0.2e1 * t399) * t365 + (t282 * t365 + (t283 - t403) * t362) * t293) * t312) * t284, t400 + 0.2e1 * (t282 * t293 * t284 + (-t284 * t417 - t293 * t416) * t295) * t295;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end