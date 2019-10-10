% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR1
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
%   Wie in S6PPRRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:18
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRRR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobiaD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:21
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
	% StartTime: 2019-10-09 21:18:21
	% EndTime: 2019-10-09 21:18:22
	% DurationCPUTime: 0.95s
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
	% StartTime: 2019-10-09 21:18:21
	% EndTime: 2019-10-09 21:18:22
	% DurationCPUTime: 1.01s
	% Computational Cost: add. (3188->68), mult. (9016->146), div. (299->12), fcn. (11669->15), ass. (0->82)
	t242 = sin(pkin(7));
	t243 = sin(pkin(6));
	t245 = cos(pkin(6));
	t241 = sin(pkin(13));
	t283 = sin(pkin(12));
	t267 = t283 * t241;
	t244 = cos(pkin(13));
	t284 = cos(pkin(12));
	t268 = t284 * t244;
	t285 = cos(pkin(7));
	t289 = (t245 * t268 - t267) * t285 - t242 * t243 * t284;
	t266 = t283 * t244;
	t269 = t284 * t241;
	t256 = t245 * t266 + t269;
	t271 = t243 * t283;
	t288 = -t242 * t271 + t256 * t285;
	t232 = t245 * t269 + t266;
	t246 = sin(qJ(3));
	t286 = cos(qJ(3));
	t218 = t232 * t286 + t289 * t246;
	t270 = t244 * t285;
	t273 = t242 * t245;
	t287 = (-t241 * t246 + t286 * t270) * t243 + t286 * t273;
	t216 = t232 * t246 - t289 * t286;
	t209 = atan2(-t216, -t287);
	t204 = sin(t209);
	t205 = cos(t209);
	t192 = -t204 * t216 - t205 * t287;
	t189 = 0.1e1 / t192;
	t233 = -t245 * t267 + t268;
	t220 = t233 * t286 - t288 * t246;
	t229 = t256 * t242 + t285 * t271;
	t240 = qJ(4) + qJ(5);
	t237 = sin(t240);
	t238 = cos(t240);
	t203 = t220 * t238 + t229 * t237;
	t199 = 0.1e1 / t203;
	t224 = 0.1e1 / t287;
	t190 = 0.1e1 / t192 ^ 2;
	t200 = 0.1e1 / t203 ^ 2;
	t225 = 0.1e1 / t287 ^ 2;
	t214 = t216 ^ 2;
	t208 = t214 * t225 + 0.1e1;
	t206 = 0.1e1 / t208;
	t211 = t218 * qJD(3);
	t228 = t246 * t273 + (t286 * t241 + t246 * t270) * t243;
	t223 = t228 * qJD(3);
	t277 = t216 * t225;
	t183 = (t211 * t224 + t223 * t277) * t206;
	t261 = t204 * t287 - t205 * t216;
	t180 = t261 * t183 - t204 * t211 + t205 * t223;
	t282 = t180 * t189 * t190;
	t202 = t220 * t237 - t229 * t238;
	t198 = t202 ^ 2;
	t195 = t198 * t200 + 0.1e1;
	t219 = t233 * t246 + t288 * t286;
	t212 = t219 * qJD(3);
	t239 = qJD(4) + qJD(5);
	t265 = t229 * t239 - t212;
	t275 = t220 * t239;
	t196 = t265 * t237 + t238 * t275;
	t278 = t200 * t202;
	t197 = -t237 * t275 + t265 * t238;
	t279 = t197 * t199 * t200;
	t281 = (t196 * t278 - t198 * t279) / t195 ^ 2;
	t280 = t190 * t219;
	t276 = t216 * t228;
	t274 = t223 * t224 * t225;
	t272 = -0.2e1 * t281;
	t259 = -t199 * t237 + t238 * t278;
	t258 = t218 * t224 + t225 * t276;
	t222 = t287 * qJD(3);
	t215 = t219 ^ 2;
	t213 = t220 * qJD(3);
	t210 = t216 * qJD(3);
	t193 = 0.1e1 / t195;
	t187 = t215 * t190 + 0.1e1;
	t184 = t258 * t206;
	t181 = t261 * t184 - t204 * t218 + t205 * t228;
	t179 = -0.2e1 * t258 / t208 ^ 2 * (t211 * t277 + t214 * t274) + (0.2e1 * t274 * t276 - t210 * t224 + (t211 * t228 + t216 * t222 + t218 * t223) * t225) * t206;
	t177 = t272 + 0.2e1 * (t193 * t196 * t200 + (-t193 * t279 - t200 * t281) * t202) * t202;
	t1 = [0, 0, t179, 0, 0, 0; 0, 0, 0.2e1 * (t181 * t280 - t189 * t220) / t187 ^ 2 * (t213 * t280 - t215 * t282) + (-t212 * t189 + (-t220 * t180 - t181 * t213) * t190 + (0.2e1 * t181 * t282 + (-(-t179 * t216 - t184 * t211 + t222 + (t184 * t287 - t218) * t183) * t205 - (t179 * t287 - t184 * t223 + t210 + (t184 * t216 - t228) * t183) * t204) * t190) * t219) / t187, 0, 0, 0; 0, 0, t259 * t219 * t272 + (t259 * t213 + ((-t199 * t239 - 0.2e1 * t202 * t279) * t238 + (t196 * t238 + (-t202 * t239 + t197) * t237) * t200) * t219) * t193, t177, t177, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:21
	% EndTime: 2019-10-09 21:18:25
	% DurationCPUTime: 3.52s
	% Computational Cost: add. (18272->124), mult. (40377->248), div. (822->12), fcn. (53065->17), ass. (0->125)
	t414 = sin(pkin(13));
	t415 = sin(pkin(12));
	t381 = t415 * t414;
	t418 = cos(pkin(13));
	t419 = cos(pkin(12));
	t387 = t419 * t418;
	t421 = cos(pkin(6));
	t348 = -t421 * t381 + t387;
	t355 = sin(qJ(3));
	t357 = cos(qJ(3));
	t383 = t415 * t418;
	t386 = t419 * t414;
	t373 = t421 * t383 + t386;
	t417 = sin(pkin(6));
	t382 = t415 * t417;
	t416 = sin(pkin(7));
	t420 = cos(pkin(7));
	t423 = t373 * t420 - t416 * t382;
	t337 = t348 * t355 + t423 * t357;
	t347 = t421 * t386 + t383;
	t372 = -t421 * t387 + t381;
	t385 = t417 * t416;
	t367 = -t372 * t420 - t419 * t385;
	t336 = t347 * t357 + t367 * t355;
	t353 = qJ(4) + qJ(5);
	t350 = sin(t353);
	t335 = -t347 * t355 + t367 * t357;
	t352 = qJD(4) + qJD(5);
	t388 = t420 * t417;
	t366 = t372 * t416 - t419 * t388;
	t365 = t335 * qJD(3) + t366 * t352;
	t351 = cos(t353);
	t397 = t351 * t352;
	t298 = t336 * t397 + t365 * t350;
	t320 = t336 * t350 - t366 * t351;
	t318 = t320 ^ 2;
	t371 = t418 * t388 + t421 * t416;
	t384 = t417 * t414;
	t344 = t371 * t355 + t357 * t384;
	t346 = -t418 * t385 + t421 * t420;
	t333 = t344 * t350 - t346 * t351;
	t330 = 0.1e1 / t333 ^ 2;
	t312 = t318 * t330 + 0.1e1;
	t310 = 0.1e1 / t312;
	t343 = -t355 * t384 + t371 * t357;
	t390 = t343 * qJD(3) + t346 * t352;
	t316 = t344 * t397 + t390 * t350;
	t329 = 0.1e1 / t333;
	t404 = t320 * t330;
	t282 = (-t298 * t329 + t316 * t404) * t310;
	t313 = atan2(-t320, t333);
	t308 = sin(t313);
	t309 = cos(t313);
	t380 = -t308 * t333 - t309 * t320;
	t278 = t380 * t282 - t298 * t308 + t309 * t316;
	t292 = -t308 * t320 + t309 * t333;
	t289 = 0.1e1 / t292;
	t290 = 0.1e1 / t292 ^ 2;
	t426 = t278 * t289 * t290;
	t338 = t348 * t357 - t355 * t423;
	t368 = t373 * t416 + t420 * t382;
	t323 = t338 * t350 - t368 * t351;
	t425 = 0.2e1 * t323 * t426;
	t376 = -t329 * t335 + t343 * t404;
	t424 = t350 * t376;
	t405 = t316 * t329 * t330;
	t422 = -0.2e1 * (t298 * t404 - t318 * t405) / t312 ^ 2;
	t324 = t338 * t351 + t368 * t350;
	t356 = cos(qJ(6));
	t354 = sin(qJ(6));
	t402 = t337 * t354;
	t307 = t324 * t356 + t402;
	t303 = 0.1e1 / t307;
	t304 = 0.1e1 / t307 ^ 2;
	t327 = t337 * qJD(3);
	t364 = t368 * t352 - t327;
	t398 = t350 * t352;
	t301 = -t338 * t398 + t364 * t351;
	t328 = t338 * qJD(3);
	t293 = t307 * qJD(6) + t301 * t354 - t328 * t356;
	t401 = t337 * t356;
	t306 = t324 * t354 - t401;
	t302 = t306 ^ 2;
	t297 = t302 * t304 + 0.1e1;
	t409 = t304 * t306;
	t396 = qJD(6) * t306;
	t294 = t301 * t356 + t328 * t354 - t396;
	t411 = t294 * t303 * t304;
	t413 = (t293 * t409 - t302 * t411) / t297 ^ 2;
	t412 = t290 * t323;
	t410 = t303 * t354;
	t408 = t306 * t356;
	t407 = t308 * t323;
	t406 = t309 * t323;
	t403 = t337 * t350;
	t319 = t323 ^ 2;
	t288 = t290 * t319 + 0.1e1;
	t300 = t338 * t397 + t364 * t350;
	t395 = 0.2e1 * (t300 * t412 - t319 * t426) / t288 ^ 2;
	t394 = -0.2e1 * t413;
	t392 = t306 * t411;
	t391 = -0.2e1 * t320 * t405;
	t389 = qJD(6) * t337 * t351 - t327;
	t378 = t304 * t408 - t410;
	t322 = t336 * t351 + t366 * t350;
	t334 = t344 * t351 + t346 * t350;
	t377 = -t322 * t329 + t334 * t404;
	t374 = qJD(6) * t338 - t328 * t351 + t337 * t398;
	t342 = t344 * qJD(3);
	t326 = t336 * qJD(3);
	t317 = -t344 * t398 + t390 * t351;
	t315 = t338 * t354 - t351 * t401;
	t314 = -t338 * t356 - t351 * t402;
	t299 = -t336 * t398 + t365 * t351;
	t295 = 0.1e1 / t297;
	t286 = 0.1e1 / t288;
	t284 = t310 * t424;
	t283 = t377 * t310;
	t280 = (-t308 * t335 + t309 * t343) * t350 + t380 * t284;
	t279 = t380 * t283 - t308 * t322 + t309 * t334;
	t276 = t377 * t422 + (t334 * t391 - t299 * t329 + (t298 * t334 + t316 * t322 + t317 * t320) * t330) * t310;
	t275 = t422 * t424 + (t376 * t397 + (t343 * t391 + t326 * t329 + (t298 * t343 + t316 * t335 - t320 * t342) * t330) * t350) * t310;
	t274 = t378 * t323 * t394 + (t378 * t300 + ((-qJD(6) * t303 - 0.2e1 * t392) * t356 + (t293 * t356 + (t294 - t396) * t354) * t304) * t323) * t295;
	t273 = (t279 * t412 - t289 * t324) * t395 + (t279 * t425 + t301 * t289 + (-t324 * t278 - t279 * t300 - (-t276 * t320 - t283 * t298 + t317 + (-t283 * t333 - t322) * t282) * t406 - (-t276 * t333 - t283 * t316 - t299 + (t283 * t320 - t334) * t282) * t407) * t290) * t286;
	t1 = [0, 0, t275, t276, t276, 0; 0, 0, (t280 * t412 + t289 * t403) * t395 + ((-t328 * t350 - t337 * t397) * t289 + t280 * t425 + (-t280 * t300 + t403 * t278 - (t343 * t397 - t275 * t320 - t284 * t298 - t342 * t350 + (-t284 * t333 - t335 * t350) * t282) * t406 - (-t335 * t397 - t275 * t333 - t284 * t316 + t326 * t350 + (t284 * t320 - t343 * t350) * t282) * t407) * t290) * t286, t273, t273, 0; 0, 0, 0.2e1 * (-t303 * t314 + t315 * t409) * t413 + (0.2e1 * t315 * t392 - t389 * t303 * t356 + t374 * t410 + (-t389 * t306 * t354 - t315 * t293 - t314 * t294 - t374 * t408) * t304) * t295, t274, t274, t394 + 0.2e1 * (t293 * t295 * t304 + (-t295 * t411 - t304 * t413) * t306) * t306;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end