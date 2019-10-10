% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR2
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
%   Wie in S6PPRRRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:20
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRRR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_jacobiaD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (151->10), mult. (481->28), div. (18->4), fcn. (595->10), ass. (0->20)
	t84 = sin(pkin(13));
	t88 = cos(pkin(13));
	t89 = cos(pkin(12));
	t85 = sin(pkin(12));
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
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:20
	% DurationCPUTime: 0.98s
	% Computational Cost: add. (2591->65), mult. (8150->143), div. (281->12), fcn. (10563->15), ass. (0->78)
	t208 = sin(pkin(7));
	t209 = sin(pkin(6));
	t211 = cos(pkin(6));
	t207 = sin(pkin(13));
	t250 = sin(pkin(12));
	t234 = t250 * t207;
	t210 = cos(pkin(13));
	t251 = cos(pkin(12));
	t235 = t251 * t210;
	t252 = cos(pkin(7));
	t256 = (t211 * t235 - t234) * t252 - t208 * t209 * t251;
	t233 = t250 * t210;
	t236 = t251 * t207;
	t224 = t211 * t233 + t236;
	t238 = t209 * t250;
	t255 = -t208 * t238 + t224 * t252;
	t202 = t211 * t236 + t233;
	t213 = sin(qJ(3));
	t253 = cos(qJ(3));
	t188 = t202 * t253 + t213 * t256;
	t237 = t210 * t252;
	t241 = t208 * t211;
	t254 = (-t207 * t213 + t253 * t237) * t209 + t253 * t241;
	t186 = t202 * t213 - t253 * t256;
	t179 = atan2(-t186, -t254);
	t174 = sin(t179);
	t175 = cos(t179);
	t162 = -t174 * t186 - t175 * t254;
	t159 = 0.1e1 / t162;
	t203 = -t211 * t234 + t235;
	t190 = t203 * t253 - t255 * t213;
	t199 = t224 * t208 + t252 * t238;
	t212 = sin(qJ(4));
	t214 = cos(qJ(4));
	t173 = t190 * t214 + t199 * t212;
	t169 = 0.1e1 / t173;
	t194 = 0.1e1 / t254;
	t160 = 0.1e1 / t162 ^ 2;
	t170 = 0.1e1 / t173 ^ 2;
	t195 = 0.1e1 / t254 ^ 2;
	t184 = t186 ^ 2;
	t178 = t184 * t195 + 0.1e1;
	t176 = 0.1e1 / t178;
	t181 = t188 * qJD(3);
	t198 = t213 * t241 + (t253 * t207 + t213 * t237) * t209;
	t192 = t198 * qJD(3);
	t244 = t186 * t195;
	t153 = (t181 * t194 + t192 * t244) * t176;
	t229 = t174 * t254 - t175 * t186;
	t150 = t229 * t153 - t174 * t181 + t175 * t192;
	t249 = t150 * t159 * t160;
	t172 = t190 * t212 - t199 * t214;
	t168 = t172 ^ 2;
	t165 = t168 * t170 + 0.1e1;
	t189 = t203 * t213 + t255 * t253;
	t182 = t189 * qJD(3);
	t166 = t173 * qJD(4) - t182 * t212;
	t245 = t170 * t172;
	t240 = qJD(4) * t172;
	t167 = -t182 * t214 - t240;
	t246 = t167 * t169 * t170;
	t248 = (t166 * t245 - t168 * t246) / t165 ^ 2;
	t247 = t160 * t189;
	t243 = t186 * t198;
	t242 = t192 * t194 * t195;
	t239 = -0.2e1 * t248;
	t227 = -t169 * t212 + t214 * t245;
	t226 = t188 * t194 + t195 * t243;
	t191 = t254 * qJD(3);
	t185 = t189 ^ 2;
	t183 = t190 * qJD(3);
	t180 = t186 * qJD(3);
	t163 = 0.1e1 / t165;
	t157 = t185 * t160 + 0.1e1;
	t154 = t226 * t176;
	t151 = t229 * t154 - t174 * t188 + t175 * t198;
	t149 = -0.2e1 * t226 / t178 ^ 2 * (t181 * t244 + t184 * t242) + (0.2e1 * t242 * t243 - t180 * t194 + (t181 * t198 + t186 * t191 + t188 * t192) * t195) * t176;
	t1 = [0, 0, t149, 0, 0, 0; 0, 0, 0.2e1 * (t151 * t247 - t159 * t190) / t157 ^ 2 * (t183 * t247 - t185 * t249) + (-t182 * t159 + (-t190 * t150 - t151 * t183) * t160 + (0.2e1 * t151 * t249 + (-(-t149 * t186 - t154 * t181 + t191 + (t154 * t254 - t188) * t153) * t175 - (t149 * t254 - t154 * t192 + t180 + (t154 * t186 - t198) * t153) * t174) * t160) * t189) / t157, 0, 0, 0; 0, 0, t227 * t189 * t239 + (t227 * t183 + ((-qJD(4) * t169 - 0.2e1 * t172 * t246) * t214 + (t166 * t214 + (t167 - t240) * t212) * t170) * t189) * t163, t239 + 0.2e1 * (t163 * t166 * t170 + (-t163 * t246 - t170 * t248) * t172) * t172, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:20
	% EndTime: 2019-10-09 21:20:22
	% DurationCPUTime: 2.48s
	% Computational Cost: add. (9293->119), mult. (27507->244), div. (559->12), fcn. (36133->17), ass. (0->119)
	t368 = sin(pkin(13));
	t369 = sin(pkin(12));
	t337 = t369 * t368;
	t372 = cos(pkin(13));
	t373 = cos(pkin(12));
	t343 = t373 * t372;
	t375 = cos(pkin(6));
	t308 = -t375 * t337 + t343;
	t312 = sin(qJ(3));
	t315 = cos(qJ(3));
	t339 = t369 * t372;
	t342 = t373 * t368;
	t329 = t375 * t339 + t342;
	t371 = sin(pkin(6));
	t338 = t369 * t371;
	t370 = sin(pkin(7));
	t374 = cos(pkin(7));
	t377 = t329 * t374 - t370 * t338;
	t292 = t308 * t312 + t377 * t315;
	t307 = t375 * t342 + t339;
	t328 = -t375 * t343 + t337;
	t341 = t371 * t370;
	t323 = -t328 * t374 - t373 * t341;
	t291 = t307 * t315 + t323 * t312;
	t311 = sin(qJ(4));
	t314 = cos(qJ(4));
	t344 = t374 * t371;
	t322 = t328 * t370 - t373 * t344;
	t282 = t291 * t314 + t322 * t311;
	t290 = -t307 * t312 + t323 * t315;
	t285 = t290 * qJD(3);
	t258 = t282 * qJD(4) + t285 * t311;
	t280 = t291 * t311 - t322 * t314;
	t278 = t280 ^ 2;
	t327 = t372 * t344 + t375 * t370;
	t340 = t371 * t368;
	t304 = t327 * t312 + t315 * t340;
	t306 = -t372 * t341 + t375 * t374;
	t297 = t304 * t311 - t306 * t314;
	t295 = 0.1e1 / t297 ^ 2;
	t272 = t278 * t295 + 0.1e1;
	t270 = 0.1e1 / t272;
	t298 = t304 * t314 + t306 * t311;
	t303 = -t312 * t340 + t327 * t315;
	t299 = t303 * qJD(3);
	t276 = t298 * qJD(4) + t299 * t311;
	t294 = 0.1e1 / t297;
	t357 = t280 * t295;
	t242 = (-t258 * t294 + t276 * t357) * t270;
	t273 = atan2(-t280, t297);
	t268 = sin(t273);
	t269 = cos(t273);
	t336 = -t268 * t297 - t269 * t280;
	t238 = t336 * t242 - t268 * t258 + t269 * t276;
	t252 = -t268 * t280 + t269 * t297;
	t249 = 0.1e1 / t252;
	t250 = 0.1e1 / t252 ^ 2;
	t380 = t238 * t249 * t250;
	t293 = t308 * t315 - t312 * t377;
	t324 = t329 * t370 + t374 * t338;
	t283 = t293 * t311 - t324 * t314;
	t379 = 0.2e1 * t283 * t380;
	t332 = -t290 * t294 + t303 * t357;
	t378 = t311 * t332;
	t358 = t276 * t294 * t295;
	t376 = -0.2e1 * (t258 * t357 - t278 * t358) / t272 ^ 2;
	t284 = t293 * t314 + t324 * t311;
	t310 = sin(qJ(5));
	t313 = cos(qJ(5));
	t267 = t284 * t313 + t292 * t310;
	t263 = 0.1e1 / t267;
	t264 = 0.1e1 / t267 ^ 2;
	t287 = t292 * qJD(3);
	t261 = -t283 * qJD(4) - t287 * t314;
	t288 = t293 * qJD(3);
	t253 = t267 * qJD(5) + t261 * t310 - t288 * t313;
	t266 = t284 * t310 - t292 * t313;
	t262 = t266 ^ 2;
	t257 = t262 * t264 + 0.1e1;
	t362 = t264 * t266;
	t351 = qJD(5) * t266;
	t254 = t261 * t313 + t288 * t310 - t351;
	t365 = t254 * t263 * t264;
	t367 = (t253 * t362 - t262 * t365) / t257 ^ 2;
	t366 = t250 * t283;
	t260 = t284 * qJD(4) - t287 * t311;
	t364 = t260 * t250;
	t363 = t263 * t310;
	t361 = t266 * t313;
	t360 = t268 * t283;
	t359 = t269 * t283;
	t356 = t292 * t311;
	t355 = t292 * t314;
	t352 = qJD(4) * t314;
	t279 = t283 ^ 2;
	t248 = t279 * t250 + 0.1e1;
	t350 = 0.2e1 * (-t279 * t380 + t283 * t364) / t248 ^ 2;
	t349 = -0.2e1 * t367;
	t347 = t266 * t365;
	t346 = -0.2e1 * t280 * t358;
	t345 = qJD(5) * t355 - t287;
	t334 = t264 * t361 - t363;
	t333 = -t282 * t294 + t298 * t357;
	t330 = qJD(4) * t356 + qJD(5) * t293 - t288 * t314;
	t300 = t304 * qJD(3);
	t286 = t291 * qJD(3);
	t277 = -t297 * qJD(4) + t299 * t314;
	t275 = t293 * t310 - t313 * t355;
	t274 = -t293 * t313 - t310 * t355;
	t259 = -t280 * qJD(4) + t285 * t314;
	t255 = 0.1e1 / t257;
	t246 = 0.1e1 / t248;
	t244 = t270 * t378;
	t243 = t333 * t270;
	t240 = (-t268 * t290 + t269 * t303) * t311 + t336 * t244;
	t239 = t336 * t243 - t268 * t282 + t269 * t298;
	t236 = t333 * t376 + (t298 * t346 - t259 * t294 + (t258 * t298 + t276 * t282 + t277 * t280) * t295) * t270;
	t235 = t376 * t378 + (t332 * t352 + (t303 * t346 + t286 * t294 + (t258 * t303 + t276 * t290 - t280 * t300) * t295) * t311) * t270;
	t1 = [0, 0, t235, t236, 0, 0; 0, 0, (t240 * t366 + t249 * t356) * t350 + ((-t288 * t311 - t292 * t352) * t249 + (-t364 + t379) * t240 + (t356 * t238 - (t303 * t352 - t235 * t280 - t244 * t258 - t300 * t311 + (-t244 * t297 - t290 * t311) * t242) * t359 - (-t290 * t352 - t235 * t297 - t244 * t276 + t286 * t311 + (t244 * t280 - t303 * t311) * t242) * t360) * t250) * t246, (t239 * t366 - t249 * t284) * t350 + (t239 * t379 + t261 * t249 + (-t284 * t238 - t239 * t260 - (-t236 * t280 - t243 * t258 + t277 + (-t243 * t297 - t282) * t242) * t359 - (-t236 * t297 - t243 * t276 - t259 + (t243 * t280 - t298) * t242) * t360) * t250) * t246, 0, 0; 0, 0, 0.2e1 * (-t263 * t274 + t275 * t362) * t367 + (0.2e1 * t275 * t347 - t345 * t263 * t313 + t330 * t363 + (-t345 * t266 * t310 - t275 * t253 - t274 * t254 - t330 * t361) * t264) * t255, t334 * t283 * t349 + (t334 * t260 + ((-qJD(5) * t263 - 0.2e1 * t347) * t313 + (t253 * t313 + (t254 - t351) * t310) * t264) * t283) * t255, t349 + 0.2e1 * (t253 * t264 * t255 + (-t255 * t365 - t264 * t367) * t266) * t266, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:20
	% EndTime: 2019-10-09 21:20:22
	% DurationCPUTime: 2.64s
	% Computational Cost: add. (10307->121), mult. (29136->245), div. (577->12), fcn. (38233->17), ass. (0->123)
	t401 = sin(pkin(13));
	t402 = sin(pkin(12));
	t369 = t402 * t401;
	t405 = cos(pkin(13));
	t406 = cos(pkin(12));
	t375 = t406 * t405;
	t408 = cos(pkin(6));
	t338 = -t408 * t369 + t375;
	t345 = sin(qJ(3));
	t347 = cos(qJ(3));
	t371 = t402 * t405;
	t374 = t406 * t401;
	t361 = t408 * t371 + t374;
	t404 = sin(pkin(6));
	t370 = t402 * t404;
	t403 = sin(pkin(7));
	t407 = cos(pkin(7));
	t410 = t361 * t407 - t403 * t370;
	t322 = t338 * t345 + t410 * t347;
	t337 = t408 * t374 + t371;
	t360 = -t408 * t375 + t369;
	t373 = t404 * t403;
	t355 = -t360 * t407 - t406 * t373;
	t321 = t337 * t347 + t355 * t345;
	t344 = sin(qJ(4));
	t346 = cos(qJ(4));
	t376 = t407 * t404;
	t354 = t360 * t403 - t406 * t376;
	t312 = t321 * t346 + t354 * t344;
	t320 = -t337 * t345 + t355 * t347;
	t316 = t320 * qJD(3);
	t294 = t312 * qJD(4) + t316 * t344;
	t310 = t321 * t344 - t354 * t346;
	t308 = t310 ^ 2;
	t359 = t405 * t376 + t408 * t403;
	t372 = t404 * t401;
	t334 = t359 * t345 + t347 * t372;
	t336 = -t405 * t373 + t408 * t407;
	t327 = t334 * t344 - t336 * t346;
	t325 = 0.1e1 / t327 ^ 2;
	t302 = t308 * t325 + 0.1e1;
	t300 = 0.1e1 / t302;
	t328 = t334 * t346 + t336 * t344;
	t333 = -t345 * t372 + t359 * t347;
	t329 = t333 * qJD(3);
	t306 = t328 * qJD(4) + t329 * t344;
	t324 = 0.1e1 / t327;
	t390 = t310 * t325;
	t272 = (-t294 * t324 + t306 * t390) * t300;
	t303 = atan2(-t310, t327);
	t298 = sin(t303);
	t299 = cos(t303);
	t368 = -t298 * t327 - t299 * t310;
	t268 = t368 * t272 - t294 * t298 + t299 * t306;
	t282 = -t298 * t310 + t299 * t327;
	t279 = 0.1e1 / t282;
	t280 = 0.1e1 / t282 ^ 2;
	t413 = t268 * t279 * t280;
	t323 = t338 * t347 - t345 * t410;
	t356 = t361 * t403 + t407 * t370;
	t313 = t323 * t344 - t356 * t346;
	t412 = 0.2e1 * t313 * t413;
	t364 = -t320 * t324 + t333 * t390;
	t411 = t344 * t364;
	t391 = t306 * t324 * t325;
	t409 = -0.2e1 * (t294 * t390 - t308 * t391) / t302 ^ 2;
	t314 = t323 * t346 + t356 * t344;
	t343 = qJ(5) + qJ(6);
	t340 = sin(t343);
	t341 = cos(t343);
	t293 = t314 * t341 + t322 * t340;
	t289 = 0.1e1 / t293;
	t290 = 0.1e1 / t293 ^ 2;
	t319 = t323 * qJD(3);
	t342 = qJD(5) + qJD(6);
	t378 = t314 * t342 - t319;
	t318 = t322 * qJD(3);
	t297 = -t313 * qJD(4) - t318 * t346;
	t379 = t322 * t342 + t297;
	t283 = t379 * t340 + t378 * t341;
	t292 = t314 * t340 - t322 * t341;
	t288 = t292 ^ 2;
	t287 = t288 * t290 + 0.1e1;
	t396 = t290 * t292;
	t284 = -t378 * t340 + t379 * t341;
	t398 = t284 * t289 * t290;
	t400 = (t283 * t396 - t288 * t398) / t287 ^ 2;
	t399 = t280 * t313;
	t397 = t289 * t340;
	t395 = t292 * t341;
	t296 = t314 * qJD(4) - t318 * t344;
	t394 = t296 * t280;
	t393 = t298 * t313;
	t392 = t299 * t313;
	t389 = t322 * t344;
	t388 = t322 * t346;
	t385 = qJD(4) * t346;
	t309 = t313 ^ 2;
	t278 = t280 * t309 + 0.1e1;
	t384 = 0.2e1 * (-t309 * t413 + t313 * t394) / t278 ^ 2;
	t383 = -0.2e1 * t400;
	t381 = t292 * t398;
	t380 = -0.2e1 * t310 * t391;
	t377 = t342 * t388 - t318;
	t366 = t290 * t395 - t397;
	t365 = -t312 * t324 + t328 * t390;
	t362 = qJD(4) * t389 - t319 * t346 + t323 * t342;
	t330 = t334 * qJD(3);
	t317 = t321 * qJD(3);
	t307 = -t327 * qJD(4) + t329 * t346;
	t305 = t323 * t340 - t341 * t388;
	t304 = -t323 * t341 - t340 * t388;
	t295 = -t310 * qJD(4) + t316 * t346;
	t285 = 0.1e1 / t287;
	t276 = 0.1e1 / t278;
	t274 = t300 * t411;
	t273 = t365 * t300;
	t270 = (-t298 * t320 + t299 * t333) * t344 + t368 * t274;
	t269 = t368 * t273 - t298 * t312 + t299 * t328;
	t266 = t365 * t409 + (t328 * t380 - t295 * t324 + (t294 * t328 + t306 * t312 + t307 * t310) * t325) * t300;
	t265 = t409 * t411 + (t364 * t385 + (t333 * t380 + t317 * t324 + (t294 * t333 + t306 * t320 - t310 * t330) * t325) * t344) * t300;
	t264 = t383 + 0.2e1 * (t283 * t290 * t285 + (-t285 * t398 - t290 * t400) * t292) * t292;
	t1 = [0, 0, t265, t266, 0, 0; 0, 0, (t270 * t399 + t279 * t389) * t384 + ((-t319 * t344 - t322 * t385) * t279 + (-t394 + t412) * t270 + (t389 * t268 - (t333 * t385 - t265 * t310 - t274 * t294 - t330 * t344 + (-t274 * t327 - t320 * t344) * t272) * t392 - (-t320 * t385 - t265 * t327 - t274 * t306 + t317 * t344 + (t274 * t310 - t333 * t344) * t272) * t393) * t280) * t276, (t269 * t399 - t279 * t314) * t384 + (t269 * t412 + t297 * t279 + (-t314 * t268 - t269 * t296 - (-t266 * t310 - t273 * t294 + t307 + (-t273 * t327 - t312) * t272) * t392 - (-t266 * t327 - t273 * t306 - t295 + (t273 * t310 - t328) * t272) * t393) * t280) * t276, 0, 0; 0, 0, 0.2e1 * (-t289 * t304 + t305 * t396) * t400 + (0.2e1 * t305 * t381 - t377 * t289 * t341 + t362 * t397 + (-t377 * t292 * t340 - t305 * t283 - t304 * t284 - t362 * t395) * t290) * t285, t366 * t313 * t383 + (t366 * t296 + ((-t289 * t342 - 0.2e1 * t381) * t341 + (t283 * t341 + (-t292 * t342 + t284) * t340) * t290) * t313) * t285, t264, t264;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end