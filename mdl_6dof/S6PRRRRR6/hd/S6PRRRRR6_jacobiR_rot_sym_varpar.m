% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:23
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobiR_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:35
	% EndTime: 2019-10-09 23:23:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(14));
	t20 = sin(pkin(6));
	t19 = sin(pkin(14));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:36
	% EndTime: 2019-10-09 23:23:36
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (38->20), mult. (118->51), div. (0->0), fcn. (169->10), ass. (0->27)
	t81 = sin(pkin(7));
	t82 = sin(pkin(6));
	t101 = t81 * t82;
	t85 = cos(pkin(6));
	t100 = t81 * t85;
	t84 = cos(pkin(7));
	t86 = sin(qJ(3));
	t99 = t84 * t86;
	t88 = cos(qJ(3));
	t98 = t84 * t88;
	t87 = sin(qJ(2));
	t97 = t85 * t87;
	t89 = cos(qJ(2));
	t96 = t85 * t89;
	t95 = t86 * t87;
	t94 = t86 * t89;
	t93 = t87 * t88;
	t92 = t88 * t89;
	t80 = sin(pkin(14));
	t83 = cos(pkin(14));
	t75 = -t80 * t87 + t83 * t96;
	t91 = t83 * t101 - t75 * t84;
	t77 = -t80 * t96 - t83 * t87;
	t90 = t80 * t101 + t77 * t84;
	t78 = -t80 * t97 + t83 * t89;
	t76 = t80 * t89 + t83 * t97;
	t1 = [0, t77 * t88 - t78 * t99, -t78 * t86 + t90 * t88, 0, 0, 0; 0, t75 * t88 - t76 * t99, -t76 * t86 - t91 * t88, 0, 0, 0; 0, (-t84 * t95 + t92) * t82, t88 * t100 + (t84 * t92 - t95) * t82, 0, 0, 0; 0, -t77 * t86 - t78 * t98, -t78 * t88 - t90 * t86, 0, 0, 0; 0, -t75 * t86 - t76 * t98, -t76 * t88 + t91 * t86, 0, 0, 0; 0, (-t84 * t93 - t94) * t82, -t86 * t100 + (-t84 * t94 - t93) * t82, 0, 0, 0; 0, t78 * t81, 0, 0, 0, 0; 0, t76 * t81, 0, 0, 0, 0; 0, t87 * t101, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:36
	% EndTime: 2019-10-09 23:23:36
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (177->51), mult. (541->118), div. (0->0), fcn. (739->14), ass. (0->55)
	t181 = sin(pkin(8));
	t182 = sin(pkin(7));
	t217 = t181 * t182;
	t183 = sin(pkin(6));
	t216 = t182 * t183;
	t185 = cos(pkin(8));
	t215 = t182 * t185;
	t187 = cos(pkin(6));
	t214 = t182 * t187;
	t186 = cos(pkin(7));
	t213 = t183 * t186;
	t188 = sin(qJ(4));
	t212 = t185 * t188;
	t191 = cos(qJ(4));
	t211 = t185 * t191;
	t189 = sin(qJ(3));
	t210 = t186 * t189;
	t192 = cos(qJ(3));
	t209 = t186 * t192;
	t190 = sin(qJ(2));
	t208 = t187 * t190;
	t193 = cos(qJ(2));
	t207 = t187 * t193;
	t206 = t189 * t190;
	t205 = t189 * t193;
	t204 = t190 * t192;
	t203 = t192 * t193;
	t184 = cos(pkin(14));
	t202 = t184 * t216;
	t201 = t190 * t216;
	t180 = sin(pkin(14));
	t174 = -t180 * t190 + t184 * t207;
	t175 = t180 * t193 + t184 * t208;
	t159 = -t175 * t189 + (t174 * t186 - t202) * t192;
	t200 = t159 * t185 + (-t174 * t182 - t184 * t213) * t181;
	t177 = -t180 * t208 + t184 * t193;
	t176 = -t180 * t207 - t184 * t190;
	t195 = t176 * t186 + t180 * t216;
	t161 = -t177 * t189 + t195 * t192;
	t199 = t161 * t185 + (-t176 * t182 + t180 * t213) * t181;
	t167 = t192 * t214 + (t186 * t203 - t206) * t183;
	t198 = t167 * t185 + (t187 * t186 - t193 * t216) * t181;
	t163 = -t174 * t189 - t175 * t209;
	t197 = t163 * t185 + t175 * t217;
	t165 = -t176 * t189 - t177 * t209;
	t196 = t165 * t185 + t177 * t217;
	t171 = (-t186 * t204 - t205) * t183;
	t194 = t171 * t185 + t181 * t201;
	t172 = (-t186 * t206 + t203) * t183;
	t168 = t189 * t214 + (t186 * t205 + t204) * t183;
	t166 = t176 * t192 - t177 * t210;
	t164 = t174 * t192 - t175 * t210;
	t162 = t177 * t192 + t195 * t189;
	t160 = t174 * t210 + t175 * t192 - t189 * t202;
	t1 = [0, t166 * t191 + t196 * t188, t161 * t191 - t162 * t212, -t162 * t188 + t199 * t191, 0, 0; 0, t164 * t191 + t197 * t188, t159 * t191 - t160 * t212, -t160 * t188 + t200 * t191, 0, 0; 0, t172 * t191 + t194 * t188, t167 * t191 - t168 * t212, -t168 * t188 + t198 * t191, 0, 0; 0, -t166 * t188 + t196 * t191, -t161 * t188 - t162 * t211, -t162 * t191 - t199 * t188, 0, 0; 0, -t164 * t188 + t197 * t191, -t159 * t188 - t160 * t211, -t160 * t191 - t200 * t188, 0, 0; 0, -t172 * t188 + t194 * t191, -t167 * t188 - t168 * t211, -t168 * t191 - t198 * t188, 0, 0; 0, -t165 * t181 + t177 * t215, t162 * t181, 0, 0, 0; 0, -t163 * t181 + t175 * t215, t160 * t181, 0, 0, 0; 0, -t171 * t181 + t185 * t201, t168 * t181, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:37
	% EndTime: 2019-10-09 23:23:37
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (391->75), mult. (1178->165), div. (0->0), fcn. (1599->16), ass. (0->80)
	t264 = sin(pkin(8));
	t265 = sin(pkin(7));
	t304 = t264 * t265;
	t271 = sin(qJ(5));
	t303 = t264 * t271;
	t275 = cos(qJ(5));
	t302 = t264 * t275;
	t266 = sin(pkin(6));
	t301 = t265 * t266;
	t268 = cos(pkin(8));
	t300 = t265 * t268;
	t270 = cos(pkin(6));
	t299 = t265 * t270;
	t269 = cos(pkin(7));
	t298 = t266 * t269;
	t272 = sin(qJ(4));
	t297 = t268 * t272;
	t276 = cos(qJ(4));
	t296 = t268 * t276;
	t273 = sin(qJ(3));
	t295 = t269 * t273;
	t277 = cos(qJ(3));
	t294 = t269 * t277;
	t274 = sin(qJ(2));
	t293 = t270 * t274;
	t278 = cos(qJ(2));
	t292 = t270 * t278;
	t291 = t273 * t274;
	t290 = t273 * t278;
	t289 = t274 * t277;
	t288 = t277 * t278;
	t267 = cos(pkin(14));
	t287 = t267 * t301;
	t286 = t274 * t301;
	t263 = sin(pkin(14));
	t257 = -t263 * t274 + t267 * t292;
	t258 = t263 * t278 + t267 * t293;
	t241 = -t258 * t273 + (t257 * t269 - t287) * t277;
	t252 = -t257 * t265 - t267 * t298;
	t285 = t241 * t268 + t252 * t264;
	t260 = -t263 * t293 + t267 * t278;
	t259 = -t263 * t292 - t267 * t274;
	t280 = t259 * t269 + t263 * t301;
	t243 = -t260 * t273 + t280 * t277;
	t253 = -t259 * t265 + t263 * t298;
	t284 = t243 * t268 + t253 * t264;
	t250 = t277 * t299 + (t269 * t288 - t291) * t266;
	t256 = t270 * t269 - t278 * t301;
	t283 = t250 * t268 + t256 * t264;
	t245 = -t257 * t273 - t258 * t294;
	t282 = t245 * t268 + t258 * t304;
	t247 = -t259 * t273 - t260 * t294;
	t281 = t247 * t268 + t260 * t304;
	t254 = (-t269 * t289 - t290) * t266;
	t279 = t254 * t268 + t264 * t286;
	t255 = (-t269 * t291 + t288) * t266;
	t251 = t273 * t299 + (t269 * t290 + t289) * t266;
	t249 = -t254 * t264 + t268 * t286;
	t248 = t259 * t277 - t260 * t295;
	t246 = t257 * t277 - t258 * t295;
	t244 = t260 * t277 + t280 * t273;
	t242 = t257 * t295 + t258 * t277 - t273 * t287;
	t240 = -t250 * t264 + t256 * t268;
	t239 = t255 * t276 + t279 * t272;
	t238 = -t247 * t264 + t260 * t300;
	t237 = -t245 * t264 + t258 * t300;
	t236 = t250 * t276 - t251 * t297;
	t235 = -t243 * t264 + t253 * t268;
	t234 = -t241 * t264 + t252 * t268;
	t233 = t251 * t276 + t283 * t272;
	t232 = -t251 * t272 + t283 * t276;
	t231 = t243 * t276 - t244 * t297;
	t230 = t241 * t276 - t242 * t297;
	t229 = t248 * t276 + t281 * t272;
	t228 = t246 * t276 + t282 * t272;
	t227 = t244 * t276 + t284 * t272;
	t226 = -t244 * t272 + t284 * t276;
	t225 = t242 * t276 + t285 * t272;
	t224 = -t242 * t272 + t285 * t276;
	t1 = [0, t229 * t275 + t238 * t271, t231 * t275 + t244 * t303, t226 * t275, -t227 * t271 + t235 * t275, 0; 0, t228 * t275 + t237 * t271, t230 * t275 + t242 * t303, t224 * t275, -t225 * t271 + t234 * t275, 0; 0, t239 * t275 + t249 * t271, t236 * t275 + t251 * t303, t232 * t275, -t233 * t271 + t240 * t275, 0; 0, -t229 * t271 + t238 * t275, -t231 * t271 + t244 * t302, -t226 * t271, -t227 * t275 - t235 * t271, 0; 0, -t228 * t271 + t237 * t275, -t230 * t271 + t242 * t302, -t224 * t271, -t225 * t275 - t234 * t271, 0; 0, -t239 * t271 + t249 * t275, -t236 * t271 + t251 * t302, -t232 * t271, -t233 * t275 - t240 * t271, 0; 0, t248 * t272 - t281 * t276, t243 * t272 + t244 * t296, t227, 0, 0; 0, t246 * t272 - t282 * t276, t241 * t272 + t242 * t296, t225, 0, 0; 0, t255 * t272 - t279 * t276, t250 * t272 + t251 * t296, t233, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:23:38
	% EndTime: 2019-10-09 23:23:38
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (813->105), mult. (2422->221), div. (0->0), fcn. (3274->18), ass. (0->102)
	t344 = sin(pkin(8));
	t345 = sin(pkin(7));
	t391 = t344 * t345;
	t352 = sin(qJ(5));
	t390 = t344 * t352;
	t357 = cos(qJ(5));
	t389 = t344 * t357;
	t346 = sin(pkin(6));
	t388 = t345 * t346;
	t348 = cos(pkin(8));
	t387 = t345 * t348;
	t350 = cos(pkin(6));
	t386 = t345 * t350;
	t349 = cos(pkin(7));
	t385 = t346 * t349;
	t353 = sin(qJ(4));
	t384 = t348 * t353;
	t358 = cos(qJ(4));
	t383 = t348 * t358;
	t354 = sin(qJ(3));
	t382 = t349 * t354;
	t359 = cos(qJ(3));
	t381 = t349 * t359;
	t355 = sin(qJ(2));
	t380 = t350 * t355;
	t360 = cos(qJ(2));
	t379 = t350 * t360;
	t351 = sin(qJ(6));
	t378 = t351 * t357;
	t377 = t354 * t355;
	t376 = t354 * t360;
	t375 = t355 * t359;
	t356 = cos(qJ(6));
	t374 = t356 * t357;
	t373 = t359 * t360;
	t347 = cos(pkin(14));
	t372 = t347 * t388;
	t371 = t355 * t388;
	t343 = sin(pkin(14));
	t337 = -t343 * t355 + t347 * t379;
	t338 = t343 * t360 + t347 * t380;
	t324 = -t337 * t354 - t338 * t381;
	t370 = t324 * t348 + t338 * t391;
	t339 = -t343 * t379 - t347 * t355;
	t340 = -t343 * t380 + t347 * t360;
	t326 = -t339 * t354 - t340 * t381;
	t369 = t326 * t348 + t340 * t391;
	t368 = -t337 * t345 - t347 * t385;
	t367 = -t339 * t345 + t343 * t385;
	t366 = t339 * t349 + t343 * t388;
	t365 = t350 * t349 - t360 * t388;
	t335 = (-t349 * t375 - t376) * t346;
	t364 = t335 * t348 + t344 * t371;
	t363 = t368 * t344;
	t362 = t367 * t344;
	t361 = t365 * t344;
	t336 = (-t349 * t377 + t373) * t346;
	t333 = t354 * t386 + (t349 * t376 + t375) * t346;
	t332 = t359 * t386 + (t349 * t373 - t377) * t346;
	t328 = -t335 * t344 + t348 * t371;
	t327 = t339 * t359 - t340 * t382;
	t325 = t337 * t359 - t338 * t382;
	t323 = t340 * t359 + t366 * t354;
	t322 = -t340 * t354 + t366 * t359;
	t321 = t337 * t382 + t338 * t359 - t354 * t372;
	t320 = -t338 * t354 + (t337 * t349 - t372) * t359;
	t319 = -t332 * t344 + t365 * t348;
	t316 = t336 * t358 + t364 * t353;
	t315 = t336 * t353 - t364 * t358;
	t314 = -t326 * t344 + t340 * t387;
	t313 = -t324 * t344 + t338 * t387;
	t312 = t332 * t358 - t333 * t384;
	t311 = t332 * t353 + t333 * t383;
	t310 = -t322 * t344 + t367 * t348;
	t309 = -t320 * t344 + t368 * t348;
	t308 = t333 * t358 + (t332 * t348 + t361) * t353;
	t307 = -t332 * t383 + t333 * t353 - t358 * t361;
	t306 = t312 * t357 + t333 * t390;
	t305 = t316 * t357 + t328 * t352;
	t304 = t322 * t358 - t323 * t384;
	t303 = t322 * t353 + t323 * t383;
	t302 = t320 * t358 - t321 * t384;
	t301 = t320 * t353 + t321 * t383;
	t300 = t327 * t358 + t369 * t353;
	t299 = t327 * t353 - t369 * t358;
	t298 = t325 * t358 + t370 * t353;
	t297 = t325 * t353 - t370 * t358;
	t296 = t323 * t358 + (t322 * t348 + t362) * t353;
	t295 = -t322 * t383 + t323 * t353 - t358 * t362;
	t294 = t321 * t358 + (t320 * t348 + t363) * t353;
	t293 = -t320 * t383 + t321 * t353 - t358 * t363;
	t292 = t308 * t357 + t319 * t352;
	t291 = -t308 * t352 + t319 * t357;
	t290 = t304 * t357 + t323 * t390;
	t289 = t302 * t357 + t321 * t390;
	t288 = t300 * t357 + t314 * t352;
	t287 = t298 * t357 + t313 * t352;
	t286 = t296 * t357 + t310 * t352;
	t285 = -t296 * t352 + t310 * t357;
	t284 = t294 * t357 + t309 * t352;
	t283 = -t294 * t352 + t309 * t357;
	t1 = [0, t288 * t356 + t299 * t351, t290 * t356 + t303 * t351, -t295 * t374 + t296 * t351, t285 * t356, -t286 * t351 + t295 * t356; 0, t287 * t356 + t297 * t351, t289 * t356 + t301 * t351, -t293 * t374 + t294 * t351, t283 * t356, -t284 * t351 + t293 * t356; 0, t305 * t356 + t315 * t351, t306 * t356 + t311 * t351, -t307 * t374 + t308 * t351, t291 * t356, -t292 * t351 + t307 * t356; 0, -t288 * t351 + t299 * t356, -t290 * t351 + t303 * t356, t295 * t378 + t296 * t356, -t285 * t351, -t286 * t356 - t295 * t351; 0, -t287 * t351 + t297 * t356, -t289 * t351 + t301 * t356, t293 * t378 + t294 * t356, -t283 * t351, -t284 * t356 - t293 * t351; 0, -t305 * t351 + t315 * t356, -t306 * t351 + t311 * t356, t307 * t378 + t308 * t356, -t291 * t351, -t292 * t356 - t307 * t351; 0, t300 * t352 - t314 * t357, t304 * t352 - t323 * t389, -t295 * t352, t286, 0; 0, t298 * t352 - t313 * t357, t302 * t352 - t321 * t389, -t293 * t352, t284, 0; 0, t316 * t352 - t328 * t357, t312 * t352 - t333 * t389, -t307 * t352, t292, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end