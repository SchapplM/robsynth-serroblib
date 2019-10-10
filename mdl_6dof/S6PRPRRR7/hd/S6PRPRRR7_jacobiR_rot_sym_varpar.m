% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:05
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobiR_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(13));
	t20 = sin(pkin(6));
	t19 = sin(pkin(13));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (18->12), mult. (56->31), div. (0->0), fcn. (81->10), ass. (0->20)
	t64 = sin(pkin(14));
	t70 = cos(pkin(7));
	t78 = t64 * t70;
	t68 = cos(pkin(14));
	t77 = t68 * t70;
	t72 = sin(qJ(2));
	t76 = t70 * t72;
	t71 = cos(pkin(6));
	t75 = t71 * t72;
	t73 = cos(qJ(2));
	t74 = t71 * t73;
	t69 = cos(pkin(13));
	t67 = sin(pkin(6));
	t66 = sin(pkin(7));
	t65 = sin(pkin(13));
	t63 = t65 * t75 - t69 * t73;
	t62 = -t65 * t74 - t69 * t72;
	t61 = -t65 * t73 - t69 * t75;
	t60 = -t65 * t72 + t69 * t74;
	t1 = [0, t62 * t68 + t63 * t78, 0, 0, 0, 0; 0, t60 * t68 + t61 * t78, 0, 0, 0, 0; 0, (-t64 * t76 + t68 * t73) * t67, 0, 0, 0, 0; 0, -t62 * t64 + t63 * t77, 0, 0, 0, 0; 0, -t60 * t64 + t61 * t77, 0, 0, 0, 0; 0, (-t64 * t73 - t68 * t76) * t67, 0, 0, 0, 0; 0, -t63 * t66, 0, 0, 0, 0; 0, -t61 * t66, 0, 0, 0, 0; 0, t67 * t72 * t66, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:03
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (118->44), mult. (365->102), div. (0->0), fcn. (498->14), ass. (0->49)
	t152 = sin(pkin(14));
	t160 = cos(pkin(7));
	t186 = t152 * t160;
	t154 = sin(pkin(8));
	t155 = sin(pkin(7));
	t185 = t154 * t155;
	t156 = sin(pkin(6));
	t184 = t155 * t156;
	t159 = cos(pkin(8));
	t183 = t155 * t159;
	t161 = cos(pkin(6));
	t182 = t155 * t161;
	t181 = t156 * t160;
	t163 = sin(qJ(2));
	t180 = t156 * t163;
	t157 = cos(pkin(14));
	t179 = t157 * t160;
	t178 = t160 * t163;
	t165 = cos(qJ(2));
	t177 = t160 * t165;
	t176 = t161 * t163;
	t175 = t161 * t165;
	t174 = t155 * t180;
	t153 = sin(pkin(13));
	t158 = cos(pkin(13));
	t147 = -t153 * t163 + t158 * t175;
	t148 = t153 * t165 + t158 * t176;
	t168 = t147 * t160 - t158 * t184;
	t173 = (-t148 * t152 + t168 * t157) * t159 + (-t147 * t155 - t158 * t181) * t154;
	t149 = -t153 * t175 - t158 * t163;
	t150 = -t153 * t176 + t158 * t165;
	t167 = t149 * t160 + t153 * t184;
	t172 = (-t150 * t152 + t167 * t157) * t159 + (-t149 * t155 + t153 * t181) * t154;
	t171 = (t157 * t182 + (-t152 * t163 + t157 * t177) * t156) * t159 + (t161 * t160 - t165 * t184) * t154;
	t136 = -t147 * t152 - t148 * t179;
	t170 = t136 * t159 + t148 * t185;
	t138 = -t149 * t152 - t150 * t179;
	t169 = t138 * t159 + t150 * t185;
	t144 = (-t152 * t165 - t157 * t178) * t156;
	t166 = t144 * t159 + t154 * t174;
	t164 = cos(qJ(4));
	t162 = sin(qJ(4));
	t145 = (-t152 * t178 + t157 * t165) * t156;
	t141 = t157 * t180 + (t156 * t177 + t182) * t152;
	t139 = t149 * t157 - t150 * t186;
	t137 = t147 * t157 - t148 * t186;
	t135 = t150 * t157 + t167 * t152;
	t133 = t148 * t157 + t168 * t152;
	t1 = [0, t139 * t164 + t169 * t162, 0, -t135 * t162 + t172 * t164, 0, 0; 0, t137 * t164 + t170 * t162, 0, -t133 * t162 + t173 * t164, 0, 0; 0, t145 * t164 + t166 * t162, 0, -t141 * t162 + t171 * t164, 0, 0; 0, -t139 * t162 + t169 * t164, 0, -t135 * t164 - t172 * t162, 0, 0; 0, -t137 * t162 + t170 * t164, 0, -t133 * t164 - t173 * t162, 0, 0; 0, -t145 * t162 + t166 * t164, 0, -t141 * t164 - t171 * t162, 0, 0; 0, -t138 * t154 + t150 * t183, 0, 0, 0, 0; 0, -t136 * t154 + t148 * t183, 0, 0, 0, 0; 0, -t144 * t154 + t159 * t174, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:03
	% EndTime: 2019-10-09 22:05:03
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (296->62), mult. (885->138), div. (0->0), fcn. (1202->16), ass. (0->72)
	t237 = sin(pkin(14));
	t245 = cos(pkin(7));
	t273 = t237 * t245;
	t239 = sin(pkin(8));
	t240 = sin(pkin(7));
	t272 = t239 * t240;
	t241 = sin(pkin(6));
	t271 = t240 * t241;
	t244 = cos(pkin(8));
	t270 = t240 * t244;
	t246 = cos(pkin(6));
	t269 = t240 * t246;
	t268 = t241 * t245;
	t249 = sin(qJ(2));
	t267 = t241 * t249;
	t242 = cos(pkin(14));
	t266 = t242 * t245;
	t265 = t245 * t249;
	t252 = cos(qJ(2));
	t264 = t245 * t252;
	t263 = t246 * t249;
	t262 = t246 * t252;
	t261 = t240 * t267;
	t238 = sin(pkin(13));
	t243 = cos(pkin(13));
	t233 = t238 * t252 + t243 * t263;
	t232 = -t238 * t249 + t243 * t262;
	t255 = t232 * t245 - t243 * t271;
	t216 = -t233 * t237 + t255 * t242;
	t227 = -t232 * t240 - t243 * t268;
	t260 = t216 * t244 + t227 * t239;
	t235 = -t238 * t263 + t243 * t252;
	t234 = -t238 * t262 - t243 * t249;
	t254 = t234 * t245 + t238 * t271;
	t218 = -t235 * t237 + t254 * t242;
	t228 = -t234 * t240 + t238 * t268;
	t259 = t218 * t244 + t228 * t239;
	t225 = t242 * t269 + (-t237 * t249 + t242 * t264) * t241;
	t231 = t246 * t245 - t252 * t271;
	t258 = t225 * t244 + t231 * t239;
	t220 = -t232 * t237 - t233 * t266;
	t257 = t220 * t244 + t233 * t272;
	t222 = -t234 * t237 - t235 * t266;
	t256 = t222 * t244 + t235 * t272;
	t229 = (-t237 * t252 - t242 * t265) * t241;
	t253 = t229 * t244 + t239 * t261;
	t251 = cos(qJ(4));
	t250 = cos(qJ(5));
	t248 = sin(qJ(4));
	t247 = sin(qJ(5));
	t230 = (-t237 * t265 + t242 * t252) * t241;
	t226 = t242 * t267 + (t241 * t264 + t269) * t237;
	t224 = -t229 * t239 + t244 * t261;
	t223 = t234 * t242 - t235 * t273;
	t221 = t232 * t242 - t233 * t273;
	t219 = t235 * t242 + t254 * t237;
	t217 = t233 * t242 + t255 * t237;
	t215 = -t225 * t239 + t231 * t244;
	t214 = t230 * t251 + t253 * t248;
	t213 = -t222 * t239 + t235 * t270;
	t212 = -t220 * t239 + t233 * t270;
	t211 = -t218 * t239 + t228 * t244;
	t210 = -t216 * t239 + t227 * t244;
	t209 = t226 * t251 + t258 * t248;
	t208 = -t226 * t248 + t258 * t251;
	t207 = t223 * t251 + t256 * t248;
	t206 = t221 * t251 + t257 * t248;
	t205 = t219 * t251 + t259 * t248;
	t204 = -t219 * t248 + t259 * t251;
	t203 = t217 * t251 + t260 * t248;
	t202 = -t217 * t248 + t260 * t251;
	t1 = [0, t207 * t250 + t213 * t247, 0, t204 * t250, -t205 * t247 + t211 * t250, 0; 0, t206 * t250 + t212 * t247, 0, t202 * t250, -t203 * t247 + t210 * t250, 0; 0, t214 * t250 + t224 * t247, 0, t208 * t250, -t209 * t247 + t215 * t250, 0; 0, -t207 * t247 + t213 * t250, 0, -t204 * t247, -t205 * t250 - t211 * t247, 0; 0, -t206 * t247 + t212 * t250, 0, -t202 * t247, -t203 * t250 - t210 * t247, 0; 0, -t214 * t247 + t224 * t250, 0, -t208 * t247, -t209 * t250 - t215 * t247, 0; 0, t223 * t248 - t256 * t251, 0, t205, 0, 0; 0, t221 * t248 - t257 * t251, 0, t203, 0, 0; 0, t230 * t248 - t253 * t251, 0, t209, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:04
	% EndTime: 2019-10-09 22:05:04
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (653->83), mult. (1935->179), div. (0->0), fcn. (2618->18), ass. (0->88)
	t307 = sin(pkin(8));
	t312 = cos(pkin(8));
	t305 = sin(pkin(14));
	t309 = sin(pkin(6));
	t310 = cos(pkin(14));
	t318 = sin(qJ(2));
	t313 = cos(pkin(7));
	t322 = cos(qJ(2));
	t345 = t313 * t322;
	t308 = sin(pkin(7));
	t314 = cos(pkin(6));
	t350 = t308 * t314;
	t326 = t310 * t350 + (-t305 * t318 + t310 * t345) * t309;
	t352 = t308 * t309;
	t333 = t314 * t313 - t322 * t352;
	t357 = t333 * t307 + t326 * t312;
	t306 = sin(pkin(13));
	t311 = cos(pkin(13));
	t344 = t314 * t318;
	t303 = -t306 * t344 + t311 * t322;
	t343 = t314 * t322;
	t302 = -t306 * t343 - t311 * t318;
	t334 = t302 * t313 + t306 * t352;
	t327 = t303 * t305 - t334 * t310;
	t349 = t309 * t313;
	t335 = -t302 * t308 + t306 * t349;
	t356 = -t335 * t307 + t327 * t312;
	t301 = t306 * t322 + t311 * t344;
	t300 = -t306 * t318 + t311 * t343;
	t336 = t300 * t313 - t311 * t352;
	t328 = t301 * t305 - t336 * t310;
	t337 = -t300 * t308 - t311 * t349;
	t355 = -t337 * t307 + t328 * t312;
	t354 = t305 * t313;
	t353 = t308 * t307;
	t351 = t308 * t312;
	t348 = t309 * t318;
	t347 = t310 * t313;
	t346 = t313 * t318;
	t315 = sin(qJ(6));
	t320 = cos(qJ(5));
	t342 = t315 * t320;
	t319 = cos(qJ(6));
	t341 = t319 * t320;
	t340 = t308 * t348;
	t288 = -t300 * t305 - t301 * t347;
	t339 = t288 * t312 + t301 * t353;
	t290 = -t302 * t305 - t303 * t347;
	t338 = t290 * t312 + t303 * t353;
	t298 = (-t305 * t322 - t310 * t346) * t309;
	t332 = t298 * t312 + t307 * t340;
	t321 = cos(qJ(4));
	t317 = sin(qJ(4));
	t316 = sin(qJ(5));
	t299 = (-t305 * t346 + t310 * t322) * t309;
	t296 = t310 * t348 + (t309 * t345 + t350) * t305;
	t292 = -t298 * t307 + t312 * t340;
	t291 = t302 * t310 - t303 * t354;
	t289 = t300 * t310 - t301 * t354;
	t287 = t303 * t310 + t334 * t305;
	t286 = t301 * t310 + t336 * t305;
	t285 = -t326 * t307 + t333 * t312;
	t282 = t299 * t321 + t332 * t317;
	t281 = t299 * t317 - t332 * t321;
	t280 = -t290 * t307 + t303 * t351;
	t279 = -t288 * t307 + t301 * t351;
	t278 = t327 * t307 + t335 * t312;
	t277 = t328 * t307 + t337 * t312;
	t276 = t296 * t321 + t357 * t317;
	t275 = t296 * t317 - t357 * t321;
	t274 = t282 * t320 + t292 * t316;
	t273 = t291 * t321 + t338 * t317;
	t272 = t291 * t317 - t338 * t321;
	t271 = t289 * t321 + t339 * t317;
	t270 = t289 * t317 - t339 * t321;
	t269 = t287 * t321 - t356 * t317;
	t268 = t287 * t317 + t356 * t321;
	t267 = t286 * t321 - t355 * t317;
	t266 = t286 * t317 + t355 * t321;
	t265 = t276 * t320 + t285 * t316;
	t264 = -t276 * t316 + t285 * t320;
	t263 = t273 * t320 + t280 * t316;
	t262 = t271 * t320 + t279 * t316;
	t261 = t269 * t320 + t278 * t316;
	t260 = -t269 * t316 + t278 * t320;
	t259 = t267 * t320 + t277 * t316;
	t258 = -t267 * t316 + t277 * t320;
	t1 = [0, t263 * t319 + t272 * t315, 0, -t268 * t341 + t269 * t315, t260 * t319, -t261 * t315 + t268 * t319; 0, t262 * t319 + t270 * t315, 0, -t266 * t341 + t267 * t315, t258 * t319, -t259 * t315 + t266 * t319; 0, t274 * t319 + t281 * t315, 0, -t275 * t341 + t276 * t315, t264 * t319, -t265 * t315 + t275 * t319; 0, -t263 * t315 + t272 * t319, 0, t268 * t342 + t269 * t319, -t260 * t315, -t261 * t319 - t268 * t315; 0, -t262 * t315 + t270 * t319, 0, t266 * t342 + t267 * t319, -t258 * t315, -t259 * t319 - t266 * t315; 0, -t274 * t315 + t281 * t319, 0, t275 * t342 + t276 * t319, -t264 * t315, -t265 * t319 - t275 * t315; 0, t273 * t316 - t280 * t320, 0, -t268 * t316, t261, 0; 0, t271 * t316 - t279 * t320, 0, -t266 * t316, t259, 0; 0, t282 * t316 - t292 * t320, 0, -t275 * t316, t265, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end