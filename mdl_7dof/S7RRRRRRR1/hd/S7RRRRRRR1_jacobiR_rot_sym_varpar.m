% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% JR_rot [9x7]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 17:10
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S7RRRRRRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiR_rot_sym_varpar: qJ has to be [7x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S7RRRRRRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiR_rot_sym_varpar: pkin has to be [4x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t7 = sin(qJ(2));
	t8 = sin(qJ(1));
	t14 = t8 * t7;
	t9 = cos(qJ(2));
	t13 = t8 * t9;
	t10 = cos(qJ(1));
	t12 = t10 * t7;
	t11 = t10 * t9;
	t1 = [-t13, -t12, 0, 0, 0, 0, 0; t11, -t14, 0, 0, 0, 0, 0; 0, t9, 0, 0, 0, 0, 0; t14, -t11, 0, 0, 0, 0, 0; -t12, -t13, 0, 0, 0, 0, 0; 0, -t7, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0, 0; t8, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (17->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t58 = sin(qJ(2));
	t59 = sin(qJ(1));
	t69 = t59 * t58;
	t60 = cos(qJ(3));
	t68 = t59 * t60;
	t57 = sin(qJ(3));
	t61 = cos(qJ(2));
	t67 = t61 * t57;
	t66 = t61 * t60;
	t62 = cos(qJ(1));
	t65 = t62 * t58;
	t64 = t62 * t60;
	t63 = t62 * t61;
	t56 = -t59 * t57 + t60 * t63;
	t55 = -t57 * t63 - t68;
	t54 = -t62 * t57 - t59 * t66;
	t53 = t59 * t67 - t64;
	t1 = [t54, -t58 * t64, t55, 0, 0, 0, 0; t56, -t58 * t68, -t53, 0, 0, 0, 0; 0, t66, -t58 * t57, 0, 0, 0, 0; t53, t57 * t65, -t56, 0, 0, 0, 0; t55, t57 * t69, t54, 0, 0, 0, 0; 0, -t67, -t58 * t60, 0, 0, 0, 0; t69, -t63, 0, 0, 0, 0, 0; -t65, -t59 * t61, 0, 0, 0, 0, 0; 0, -t58, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (37->22), mult. (118->39), div. (0->0), fcn. (180->8), ass. (0->30)
	t90 = sin(qJ(4));
	t92 = sin(qJ(2));
	t112 = t92 * t90;
	t94 = cos(qJ(4));
	t111 = t92 * t94;
	t95 = cos(qJ(3));
	t110 = t92 * t95;
	t97 = cos(qJ(1));
	t109 = t92 * t97;
	t91 = sin(qJ(3));
	t93 = sin(qJ(1));
	t108 = t93 * t91;
	t107 = t93 * t95;
	t96 = cos(qJ(2));
	t106 = t96 * t90;
	t105 = t96 * t91;
	t104 = t96 * t94;
	t103 = t97 * t91;
	t102 = t97 * t95;
	t86 = t96 * t107 + t103;
	t101 = t93 * t111 - t86 * t90;
	t100 = -t93 * t112 - t86 * t94;
	t99 = -t94 * t110 + t106;
	t98 = t90 * t110 + t104;
	t88 = t96 * t102 - t108;
	t87 = -t96 * t103 - t107;
	t85 = t93 * t105 - t102;
	t84 = t90 * t109 + t88 * t94;
	t83 = t94 * t109 - t88 * t90;
	t1 = [t100, t99 * t97, t87 * t94, t83, 0, 0, 0; t84, t99 * t93, -t85 * t94, t101, 0, 0, 0; 0, t95 * t104 + t112, -t91 * t111, -t98, 0, 0, 0; -t101, t98 * t97, -t87 * t90, -t84, 0, 0, 0; t83, t98 * t93, t85 * t90, t100, 0, 0, 0; 0, -t95 * t106 + t111, t91 * t112, t99, 0, 0, 0; t85, t92 * t103, -t88, 0, 0, 0, 0; t87, t92 * t108, -t86, 0, 0, 0, 0; 0, -t105, -t110, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:05
	% EndTime: 2019-10-10 17:10:05
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (90->37), mult. (280->82), div. (0->0), fcn. (405->10), ass. (0->43)
	t166 = cos(qJ(3));
	t161 = sin(qJ(3));
	t168 = cos(qJ(1));
	t172 = t168 * t161;
	t163 = sin(qJ(1));
	t167 = cos(qJ(2));
	t176 = t163 * t167;
	t151 = t166 * t176 + t172;
	t165 = cos(qJ(4));
	t160 = sin(qJ(4));
	t162 = sin(qJ(2));
	t180 = t162 * t160;
	t142 = t151 * t165 + t163 * t180;
	t171 = t168 * t166;
	t150 = t161 * t176 - t171;
	t159 = sin(qJ(5));
	t164 = cos(qJ(5));
	t187 = t142 * t159 + t150 * t164;
	t186 = -t142 * t164 + t150 * t159;
	t183 = t159 * t165;
	t182 = t161 * t162;
	t181 = t161 * t167;
	t179 = t162 * t165;
	t178 = t162 * t166;
	t177 = t162 * t168;
	t175 = t164 * t165;
	t174 = t167 * t160;
	t173 = t167 * t165;
	t170 = t159 * t182;
	t169 = t164 * t182;
	t141 = -t151 * t160 + t163 * t179;
	t149 = t165 * t178 - t174;
	t148 = -t160 * t178 - t173;
	t154 = -t163 * t161 + t167 * t171;
	t153 = -t163 * t166 - t167 * t172;
	t152 = t166 * t173 + t180;
	t147 = t149 * t168;
	t146 = t149 * t163;
	t145 = t154 * t165 + t160 * t177;
	t144 = t154 * t160 - t165 * t177;
	t140 = t145 * t164 + t153 * t159;
	t139 = -t145 * t159 + t153 * t164;
	t1 = [t186, -t147 * t164 + t168 * t170, t153 * t175 - t154 * t159, -t144 * t164, t139, 0, 0; t140, -t146 * t164 + t163 * t170, -t150 * t175 - t151 * t159, t141 * t164, -t187, 0, 0; 0, t152 * t164 - t159 * t181, (-t159 * t166 - t161 * t175) * t162, t148 * t164, -t149 * t159 - t169, 0, 0; t187, t147 * t159 + t168 * t169, -t153 * t183 - t154 * t164, t144 * t159, -t140, 0, 0; t139, t146 * t159 + t163 * t169, t150 * t183 - t151 * t164, -t141 * t159, t186, 0, 0; 0, -t152 * t159 - t164 * t181, (t161 * t183 - t164 * t166) * t162, -t148 * t159, -t149 * t164 + t170, 0, 0; t141, t148 * t168, t153 * t160, t145, 0, 0, 0; t144, t148 * t163, -t150 * t160, t142, 0, 0, 0; 0, t166 * t174 - t179, -t161 * t180, t149, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:05
	% EndTime: 2019-10-10 17:10:06
	% DurationCPUTime: 0.55s
	% Computational Cost: add. (197->62), mult. (593->135), div. (0->0), fcn. (835->12), ass. (0->65)
	t237 = cos(qJ(3));
	t231 = sin(qJ(3));
	t239 = cos(qJ(1));
	t244 = t239 * t231;
	t233 = sin(qJ(1));
	t238 = cos(qJ(2));
	t249 = t233 * t238;
	t219 = t237 * t249 + t244;
	t236 = cos(qJ(4));
	t230 = sin(qJ(4));
	t232 = sin(qJ(2));
	t253 = t232 * t230;
	t204 = t219 * t236 + t233 * t253;
	t243 = t239 * t237;
	t218 = t231 * t249 - t243;
	t229 = sin(qJ(5));
	t235 = cos(qJ(5));
	t193 = t204 * t235 - t218 * t229;
	t252 = t232 * t236;
	t203 = t219 * t230 - t233 * t252;
	t228 = sin(qJ(6));
	t234 = cos(qJ(6));
	t265 = t193 * t228 - t203 * t234;
	t264 = -t193 * t234 - t203 * t228;
	t192 = -t204 * t229 - t218 * t235;
	t259 = t228 * t230;
	t258 = t228 * t235;
	t257 = t229 * t236;
	t256 = t230 * t234;
	t255 = t231 * t232;
	t254 = t231 * t238;
	t251 = t232 * t237;
	t250 = t232 * t239;
	t248 = t234 * t235;
	t247 = t235 * t236;
	t246 = t238 * t230;
	t245 = t238 * t236;
	t242 = t231 * t253;
	t241 = t229 * t255;
	t240 = t235 * t255;
	t217 = t236 * t251 - t246;
	t216 = t230 * t251 + t245;
	t223 = -t233 * t231 + t238 * t243;
	t222 = -t233 * t237 - t238 * t244;
	t221 = t237 * t245 + t253;
	t220 = t237 * t246 - t252;
	t214 = t217 * t239;
	t213 = t216 * t239;
	t212 = t217 * t233;
	t211 = t216 * t233;
	t210 = (-t229 * t237 - t231 * t247) * t232;
	t209 = t223 * t236 + t230 * t250;
	t208 = t223 * t230 - t236 * t250;
	t207 = t221 * t235 - t229 * t254;
	t202 = t217 * t235 - t241;
	t201 = -t217 * t229 - t240;
	t200 = -t214 * t235 + t239 * t241;
	t199 = -t212 * t235 + t233 * t241;
	t198 = t222 * t247 - t223 * t229;
	t197 = -t218 * t247 - t219 * t229;
	t196 = t209 * t235 + t222 * t229;
	t195 = t209 * t229 - t222 * t235;
	t191 = t196 * t234 + t208 * t228;
	t190 = -t196 * t228 + t208 * t234;
	t1 = [t264, t200 * t234 - t213 * t228, t198 * t234 + t222 * t259, -t208 * t248 + t209 * t228, -t195 * t234, t190, 0; t191, t199 * t234 - t211 * t228, t197 * t234 - t218 * t259, -t203 * t248 + t204 * t228, t192 * t234, -t265, 0; 0, t207 * t234 + t220 * t228, t210 * t234 - t228 * t242, -t216 * t248 + t217 * t228, t201 * t234, -t202 * t228 + t216 * t234, 0; t265, -t200 * t228 - t213 * t234, -t198 * t228 + t222 * t256, t208 * t258 + t209 * t234, t195 * t228, -t191, 0; t190, -t199 * t228 - t211 * t234, -t197 * t228 - t218 * t256, t203 * t258 + t204 * t234, -t192 * t228, t264, 0; 0, -t207 * t228 + t220 * t234, -t210 * t228 - t234 * t242, t216 * t258 + t217 * t234, -t201 * t228, -t202 * t234 - t216 * t228, 0; t192, -t214 * t229 - t239 * t240, t222 * t257 + t223 * t235, -t208 * t229, t196, 0, 0; t195, -t212 * t229 - t233 * t240, -t218 * t257 + t219 * t235, -t203 * t229, t193, 0, 0; 0, t221 * t229 + t235 * t254, (-t231 * t257 + t235 * t237) * t232, -t216 * t229, t202, 0, 0;];
	JR_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiR_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:07
	% EndTime: 2019-10-10 17:10:08
	% DurationCPUTime: 1.15s
	% Computational Cost: add. (405->92), mult. (1175->199), div. (0->0), fcn. (1631->14), ass. (0->92)
	t343 = cos(qJ(3));
	t336 = sin(qJ(3));
	t345 = cos(qJ(1));
	t350 = t345 * t336;
	t338 = sin(qJ(1));
	t344 = cos(qJ(2));
	t356 = t338 * t344;
	t323 = t343 * t356 + t350;
	t342 = cos(qJ(4));
	t335 = sin(qJ(4));
	t337 = sin(qJ(2));
	t360 = t337 * t335;
	t305 = t323 * t342 + t338 * t360;
	t349 = t345 * t343;
	t322 = t336 * t356 - t349;
	t334 = sin(qJ(5));
	t341 = cos(qJ(5));
	t286 = t305 * t341 - t322 * t334;
	t359 = t337 * t342;
	t304 = t323 * t335 - t338 * t359;
	t333 = sin(qJ(6));
	t340 = cos(qJ(6));
	t273 = t286 * t340 + t304 * t333;
	t285 = t305 * t334 + t322 * t341;
	t332 = sin(qJ(7));
	t339 = cos(qJ(7));
	t378 = t273 * t332 + t285 * t339;
	t377 = -t273 * t339 + t285 * t332;
	t374 = t286 * t333 - t304 * t340;
	t369 = t332 * t334;
	t368 = t332 * t340;
	t367 = t333 * t335;
	t366 = t333 * t341;
	t365 = t334 * t339;
	t364 = t334 * t342;
	t363 = t335 * t340;
	t362 = t336 * t337;
	t361 = t336 * t344;
	t358 = t337 * t343;
	t357 = t337 * t345;
	t355 = t339 * t340;
	t354 = t340 * t341;
	t353 = t341 * t342;
	t352 = t344 * t335;
	t351 = t344 * t342;
	t348 = t336 * t360;
	t347 = t334 * t362;
	t346 = t341 * t362;
	t321 = t342 * t358 - t352;
	t320 = t335 * t358 + t351;
	t327 = -t338 * t336 + t344 * t349;
	t326 = -t338 * t343 - t344 * t350;
	t325 = t343 * t351 + t360;
	t324 = t343 * t352 - t359;
	t317 = t321 * t345;
	t316 = t320 * t345;
	t315 = t321 * t338;
	t314 = t320 * t338;
	t313 = (-t334 * t343 - t336 * t353) * t337;
	t312 = (-t336 * t364 + t341 * t343) * t337;
	t311 = t327 * t342 + t335 * t357;
	t310 = t327 * t335 - t342 * t357;
	t309 = t325 * t341 - t334 * t361;
	t308 = t325 * t334 + t341 * t361;
	t303 = t321 * t341 - t347;
	t302 = t321 * t334 + t346;
	t301 = -t317 * t341 + t345 * t347;
	t300 = -t317 * t334 - t345 * t346;
	t299 = -t315 * t341 + t338 * t347;
	t298 = -t315 * t334 - t338 * t346;
	t297 = t313 * t340 - t333 * t348;
	t296 = t326 * t353 - t327 * t334;
	t295 = t326 * t364 + t327 * t341;
	t294 = -t322 * t353 - t323 * t334;
	t293 = -t322 * t364 + t323 * t341;
	t292 = -t320 * t354 + t321 * t333;
	t291 = t311 * t341 + t326 * t334;
	t290 = t311 * t334 - t326 * t341;
	t289 = t309 * t340 + t324 * t333;
	t284 = t303 * t340 + t320 * t333;
	t283 = -t303 * t333 + t320 * t340;
	t282 = t301 * t340 - t316 * t333;
	t281 = t299 * t340 - t314 * t333;
	t280 = t296 * t340 + t326 * t367;
	t279 = t294 * t340 - t322 * t367;
	t278 = -t310 * t354 + t311 * t333;
	t277 = -t304 * t354 + t305 * t333;
	t276 = t291 * t340 + t310 * t333;
	t275 = -t291 * t333 + t310 * t340;
	t271 = t276 * t339 - t290 * t332;
	t270 = -t276 * t332 - t290 * t339;
	t1 = [t377, t282 * t339 - t300 * t332, t280 * t339 - t295 * t332, t278 * t339 + t310 * t369, -t290 * t355 - t291 * t332, t275 * t339, t270; t271, t281 * t339 - t298 * t332, t279 * t339 - t293 * t332, t277 * t339 + t304 * t369, -t285 * t355 - t286 * t332, -t374 * t339, -t378; 0, t289 * t339 - t308 * t332, t297 * t339 - t312 * t332, t292 * t339 + t320 * t369, -t302 * t355 - t303 * t332, t283 * t339, -t284 * t332 - t302 * t339; t378, -t282 * t332 - t300 * t339, -t280 * t332 - t295 * t339, -t278 * t332 + t310 * t365, t290 * t368 - t291 * t339, -t275 * t332, -t271; t270, -t281 * t332 - t298 * t339, -t279 * t332 - t293 * t339, -t277 * t332 + t304 * t365, t285 * t368 - t286 * t339, t374 * t332, t377; 0, -t289 * t332 - t308 * t339, -t297 * t332 - t312 * t339, -t292 * t332 + t320 * t365, t302 * t368 - t303 * t339, -t283 * t332, -t284 * t339 + t302 * t332; t374, -t301 * t333 - t316 * t340, -t296 * t333 + t326 * t363, t310 * t366 + t311 * t340, t290 * t333, -t276, 0; t275, -t299 * t333 - t314 * t340, -t294 * t333 - t322 * t363, t304 * t366 + t305 * t340, t285 * t333, -t273, 0; 0, -t309 * t333 + t324 * t340, -t313 * t333 - t340 * t348, t320 * t366 + t321 * t340, t302 * t333, -t284, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,7);
end