% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:16
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR12_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR12_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobiR_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
	t38 = sin(pkin(14));
	t42 = sin(qJ(1));
	t47 = t42 * t38;
	t40 = cos(pkin(14));
	t46 = t42 * t40;
	t43 = cos(qJ(1));
	t45 = t43 * t38;
	t44 = t43 * t40;
	t41 = cos(pkin(6));
	t39 = sin(pkin(6));
	t1 = [-t41 * t45 - t46, 0, 0, 0, 0, 0; -t41 * t47 + t44, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t41 * t44 + t47, 0, 0, 0, 0, 0; -t41 * t46 - t45, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t43 * t39, 0, 0, 0, 0, 0; t42 * t39, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:46
	% EndTime: 2019-10-10 09:15:46
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (40->17), mult. (122->36), div. (0->0), fcn. (174->10), ass. (0->27)
	t103 = cos(qJ(1));
	t96 = sin(pkin(6));
	t111 = t103 * t96;
	t101 = sin(qJ(1));
	t99 = cos(pkin(6));
	t110 = t103 * t99;
	t94 = sin(pkin(14));
	t97 = cos(pkin(14));
	t88 = t101 * t94 - t97 * t110;
	t95 = sin(pkin(7));
	t98 = cos(pkin(7));
	t117 = t95 * t111 + t88 * t98;
	t100 = sin(qJ(3));
	t102 = cos(qJ(3));
	t89 = t101 * t97 + t94 * t110;
	t116 = t89 * t100 + t117 * t102;
	t114 = t94 * t96;
	t113 = t101 * t96;
	t112 = t101 * t99;
	t107 = t96 * t97 * t98 + t95 * t99;
	t90 = -t103 * t94 - t97 * t112;
	t106 = t95 * t113 + t90 * t98;
	t104 = t117 * t100 - t89 * t102;
	t91 = t103 * t97 - t94 * t112;
	t87 = t106 * t100 + t91 * t102;
	t86 = -t91 * t100 + t106 * t102;
	t1 = [t104, 0, t86, 0, 0, 0; t87, 0, -t116, 0, 0, 0; 0, 0, -t100 * t114 + t107 * t102, 0, 0, 0; t116, 0, -t87, 0, 0, 0; t86, 0, t104, 0, 0, 0; 0, 0, -t107 * t100 - t102 * t114, 0, 0, 0; t98 * t111 - t88 * t95, 0, 0, 0, 0, 0; t98 * t113 - t90 * t95, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:47
	% EndTime: 2019-10-10 09:15:47
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (183->38), mult. (540->80), div. (0->0), fcn. (741->14), ass. (0->46)
	t184 = cos(pkin(6));
	t181 = cos(pkin(14));
	t190 = cos(qJ(1));
	t198 = t190 * t181;
	t177 = sin(pkin(14));
	t187 = sin(qJ(1));
	t201 = t187 * t177;
	t171 = -t184 * t198 + t201;
	t199 = t190 * t177;
	t200 = t187 * t181;
	t172 = t184 * t199 + t200;
	t186 = sin(qJ(3));
	t189 = cos(qJ(3));
	t179 = sin(pkin(7));
	t180 = sin(pkin(6));
	t205 = t180 * t190;
	t197 = t179 * t205;
	t183 = cos(pkin(7));
	t202 = t183 * t186;
	t162 = t171 * t202 - t172 * t189 + t186 * t197;
	t185 = sin(qJ(4));
	t188 = cos(qJ(4));
	t161 = (t171 * t183 + t197) * t189 + t172 * t186;
	t167 = -t171 * t179 + t183 * t205;
	t178 = sin(pkin(8));
	t182 = cos(pkin(8));
	t195 = t161 * t182 + t167 * t178;
	t213 = t162 * t188 + t195 * t185;
	t212 = -t162 * t185 + t195 * t188;
	t207 = t179 * t184;
	t206 = t180 * t187;
	t204 = t182 * t185;
	t203 = t182 * t188;
	t174 = -t184 * t201 + t198;
	t173 = -t184 * t200 - t199;
	t191 = t173 * t183 + t179 * t206;
	t163 = -t174 * t186 + t191 * t189;
	t169 = -t173 * t179 + t183 * t206;
	t194 = t163 * t182 + t169 * t178;
	t165 = t189 * t207 + (t181 * t183 * t189 - t177 * t186) * t180;
	t193 = t165 * t182 + (-t180 * t181 * t179 + t184 * t183) * t178;
	t166 = t186 * t207 + (t177 * t189 + t181 * t202) * t180;
	t164 = t174 * t189 + t191 * t186;
	t158 = t164 * t188 + t194 * t185;
	t157 = -t164 * t185 + t194 * t188;
	t1 = [t213, 0, t163 * t188 - t164 * t204, t157, 0, 0; t158, 0, -t161 * t188 + t162 * t204, -t212, 0, 0; 0, 0, t165 * t188 - t166 * t204, -t166 * t185 + t193 * t188, 0, 0; t212, 0, -t163 * t185 - t164 * t203, -t158, 0, 0; t157, 0, t161 * t185 + t162 * t203, t213, 0, 0; 0, 0, -t165 * t185 - t166 * t203, -t166 * t188 - t193 * t185, 0, 0; -t161 * t178 + t167 * t182, 0, t164 * t178, 0, 0, 0; -t163 * t178 + t169 * t182, 0, -t162 * t178, 0, 0, 0; 0, 0, t166 * t178, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:48
	% EndTime: 2019-10-10 09:15:48
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (399->55), mult. (1173->113), div. (0->0), fcn. (1599->16), ass. (0->62)
	t259 = cos(pkin(14));
	t262 = cos(pkin(6));
	t270 = cos(qJ(1));
	t278 = t262 * t270;
	t255 = sin(pkin(14));
	t266 = sin(qJ(1));
	t287 = t255 * t266;
	t249 = -t259 * t278 + t287;
	t277 = t266 * t259;
	t250 = t255 * t278 + t277;
	t265 = sin(qJ(3));
	t269 = cos(qJ(3));
	t257 = sin(pkin(7));
	t258 = sin(pkin(6));
	t282 = t258 * t270;
	t276 = t257 * t282;
	t261 = cos(pkin(7));
	t279 = t261 * t265;
	t239 = t249 * t279 - t250 * t269 + t265 * t276;
	t264 = sin(qJ(4));
	t268 = cos(qJ(4));
	t238 = (t249 * t261 + t276) * t269 + t250 * t265;
	t245 = -t249 * t257 + t261 * t282;
	t256 = sin(pkin(8));
	t260 = cos(pkin(8));
	t274 = t238 * t260 + t245 * t256;
	t223 = t239 * t268 + t274 * t264;
	t230 = t238 * t256 - t245 * t260;
	t263 = sin(qJ(5));
	t267 = cos(qJ(5));
	t297 = t223 * t263 + t230 * t267;
	t296 = t223 * t267 - t230 * t263;
	t221 = t239 * t264 - t274 * t268;
	t251 = -t255 * t270 - t262 * t277;
	t283 = t258 * t266;
	t247 = -t251 * t257 + t261 * t283;
	t289 = t247 * t256;
	t286 = t256 * t263;
	t285 = t256 * t267;
	t284 = t257 * t262;
	t281 = t260 * t264;
	t280 = t260 * t268;
	t243 = t269 * t284 + (t259 * t261 * t269 - t255 * t265) * t258;
	t248 = -t257 * t258 * t259 + t261 * t262;
	t273 = t243 * t260 + t248 * t256;
	t271 = t251 * t261 + t257 * t283;
	t252 = t259 * t270 - t262 * t287;
	t244 = t265 * t284 + (t255 * t269 + t259 * t279) * t258;
	t241 = t252 * t269 + t271 * t265;
	t240 = -t252 * t265 + t271 * t269;
	t235 = -t243 * t256 + t248 * t260;
	t233 = t243 * t268 - t244 * t281;
	t232 = -t240 * t256 + t247 * t260;
	t229 = t244 * t268 + t273 * t264;
	t228 = -t244 * t264 + t273 * t268;
	t227 = t240 * t268 - t241 * t281;
	t226 = -t238 * t268 + t239 * t281;
	t225 = t241 * t268 + (t240 * t260 + t289) * t264;
	t224 = -t240 * t280 + t241 * t264 - t268 * t289;
	t220 = t225 * t267 + t232 * t263;
	t219 = -t225 * t263 + t232 * t267;
	t1 = [t296, 0, t227 * t267 + t241 * t286, -t224 * t267, t219, 0; t220, 0, t226 * t267 - t239 * t286, t221 * t267, t297, 0; 0, 0, t233 * t267 + t244 * t286, t228 * t267, -t229 * t263 + t235 * t267, 0; -t297, 0, -t227 * t263 + t241 * t285, t224 * t263, -t220, 0; t219, 0, -t226 * t263 - t239 * t285, -t221 * t263, t296, 0; 0, 0, -t233 * t263 + t244 * t285, -t228 * t263, -t229 * t267 - t235 * t263, 0; t221, 0, t240 * t264 + t241 * t280, t225, 0, 0; t224, 0, -t238 * t264 - t239 * t280, -t223, 0, 0; 0, 0, t243 * t264 + t244 * t280, t229, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:15:49
	% EndTime: 2019-10-10 09:15:50
	% DurationCPUTime: 1.00s
	% Computational Cost: add. (824->80), mult. (2410->158), div. (0->0), fcn. (3270->18), ass. (0->82)
	t344 = cos(pkin(6));
	t341 = cos(pkin(14));
	t354 = cos(qJ(1));
	t365 = t354 * t341;
	t337 = sin(pkin(14));
	t349 = sin(qJ(1));
	t369 = t349 * t337;
	t331 = -t344 * t365 + t369;
	t366 = t354 * t337;
	t368 = t349 * t341;
	t332 = t344 * t366 + t368;
	t348 = sin(qJ(3));
	t353 = cos(qJ(3));
	t339 = sin(pkin(7));
	t340 = sin(pkin(6));
	t374 = t340 * t354;
	t364 = t339 * t374;
	t343 = cos(pkin(7));
	t371 = t343 * t348;
	t321 = t331 * t371 - t332 * t353 + t348 * t364;
	t347 = sin(qJ(4));
	t352 = cos(qJ(4));
	t320 = (t331 * t343 + t364) * t353 + t332 * t348;
	t342 = cos(pkin(8));
	t329 = -t331 * t339 + t343 * t374;
	t338 = sin(pkin(8));
	t381 = t329 * t338;
	t363 = t320 * t342 + t381;
	t300 = t321 * t352 + t363 * t347;
	t311 = t320 * t338 - t329 * t342;
	t346 = sin(qJ(5));
	t351 = cos(qJ(5));
	t290 = t300 * t351 - t311 * t346;
	t345 = sin(qJ(6));
	t388 = t290 * t345;
	t350 = cos(qJ(6));
	t387 = t290 * t350;
	t288 = t300 * t346 + t311 * t351;
	t383 = t321 * t347;
	t378 = t338 * t346;
	t377 = t338 * t351;
	t376 = t339 * t344;
	t375 = t340 * t349;
	t373 = t342 * t347;
	t372 = t342 * t352;
	t370 = t345 * t351;
	t367 = t350 * t351;
	t333 = -t344 * t368 - t366;
	t361 = -t333 * t339 + t343 * t375;
	t360 = t333 * t343 + t339 * t375;
	t359 = -t340 * t341 * t339 + t344 * t343;
	t357 = t361 * t338;
	t356 = t359 * t338;
	t334 = -t344 * t369 + t365;
	t322 = -t334 * t348 + t360 * t353;
	t355 = -t322 * t338 + t361 * t342;
	t328 = t348 * t376 + (t337 * t353 + t341 * t371) * t340;
	t327 = t353 * t376 + (t341 * t343 * t353 - t337 * t348) * t340;
	t323 = t334 * t353 + t360 * t348;
	t317 = -t327 * t338 + t359 * t342;
	t314 = t327 * t352 - t328 * t373;
	t313 = t327 * t347 + t328 * t372;
	t309 = t328 * t352 + (t327 * t342 + t356) * t347;
	t308 = -t327 * t372 + t328 * t347 - t352 * t356;
	t307 = t314 * t351 + t328 * t378;
	t306 = t322 * t352 - t323 * t373;
	t305 = t322 * t347 + t323 * t372;
	t304 = -t320 * t352 + t321 * t373;
	t303 = -t320 * t347 - t321 * t372;
	t302 = t323 * t352 + (t322 * t342 + t357) * t347;
	t301 = -t322 * t372 + t323 * t347 - t352 * t357;
	t299 = -t363 * t352 + t383;
	t297 = t320 * t372 + t352 * t381 - t383;
	t296 = t309 * t351 + t317 * t346;
	t295 = -t309 * t346 + t317 * t351;
	t294 = t306 * t351 + t323 * t378;
	t293 = t304 * t351 - t321 * t378;
	t292 = t302 * t351 + t355 * t346;
	t291 = t302 * t346 - t355 * t351;
	t287 = t292 * t350 + t301 * t345;
	t286 = -t292 * t345 + t301 * t350;
	t1 = [t299 * t345 + t387, 0, t294 * t350 + t305 * t345, -t301 * t367 + t302 * t345, -t291 * t350, t286; t287, 0, t293 * t350 + t303 * t345, -t297 * t367 - t300 * t345, t288 * t350, t297 * t350 + t388; 0, 0, t307 * t350 + t313 * t345, -t308 * t367 + t309 * t345, t295 * t350, -t296 * t345 + t308 * t350; t299 * t350 - t388, 0, -t294 * t345 + t305 * t350, t301 * t370 + t302 * t350, t291 * t345, -t287; t286, 0, -t293 * t345 + t303 * t350, t297 * t370 - t300 * t350, -t288 * t345, -t297 * t345 + t387; 0, 0, -t307 * t345 + t313 * t350, t308 * t370 + t309 * t350, -t295 * t345, -t296 * t350 - t308 * t345; t288, 0, t306 * t346 - t323 * t377, -t301 * t346, t292, 0; t291, 0, t304 * t346 + t321 * t377, -t297 * t346, -t290, 0; 0, 0, t314 * t346 - t328 * t377, -t308 * t346, t296, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end