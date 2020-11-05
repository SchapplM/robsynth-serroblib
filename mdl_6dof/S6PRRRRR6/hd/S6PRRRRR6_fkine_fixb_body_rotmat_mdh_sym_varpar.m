% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:22
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:22:29
	% EndTime: 2020-11-04 21:22:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:22:29
	% EndTime: 2020-11-04 21:22:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t132 = cos(pkin(14));
	t131 = sin(pkin(14));
	t1 = [t132, -t131, 0, 0; t131, t132, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:22:29
	% EndTime: 2020-11-04 21:22:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t133 = sin(pkin(14));
	t134 = sin(pkin(6));
	t142 = t133 * t134;
	t135 = cos(pkin(14));
	t141 = t135 * t134;
	t136 = cos(pkin(6));
	t137 = sin(qJ(2));
	t140 = t136 * t137;
	t138 = cos(qJ(2));
	t139 = t136 * t138;
	t1 = [-t133 * t140 + t135 * t138, -t133 * t139 - t135 * t137, t142, t135 * pkin(1) + pkin(9) * t142 + 0; t133 * t138 + t135 * t140, -t133 * t137 + t135 * t139, -t141, t133 * pkin(1) - pkin(9) * t141 + 0; t134 * t137, t134 * t138, t136, t136 * pkin(9) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:22:29
	% EndTime: 2020-11-04 21:22:30
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->35), mult. (112->66), div. (0->0), fcn. (146->10), ass. (0->30)
	t149 = sin(pkin(7));
	t171 = pkin(10) * t149;
	t148 = sin(pkin(14));
	t170 = t148 * pkin(2);
	t151 = cos(pkin(14));
	t169 = t151 * pkin(2);
	t153 = cos(pkin(6));
	t168 = t149 * t153;
	t157 = cos(qJ(2));
	t167 = t149 * t157;
	t152 = cos(pkin(7));
	t147 = t152 * pkin(10) + pkin(9);
	t150 = sin(pkin(6));
	t166 = t150 * t147;
	t165 = t150 * t152;
	t155 = sin(qJ(2));
	t164 = t152 * t155;
	t163 = t152 * t157;
	t162 = t153 * t155;
	t161 = t153 * t157;
	t160 = t148 * t171;
	t159 = t151 * t171;
	t158 = -t149 * t150 + t152 * t161;
	t156 = cos(qJ(3));
	t154 = sin(qJ(3));
	t146 = t148 * t157 + t151 * t162;
	t145 = t148 * t162 - t151 * t157;
	t144 = -t148 * t164 + t158 * t151;
	t143 = -t158 * t148 - t151 * t164;
	t1 = [t143 * t154 - t156 * t145, t143 * t156 + t154 * t145, (t148 * t161 + t151 * t155) * t149 + t148 * t165, (t153 * t160 + t169) * t157 + (-t153 * t170 + t159) * t155 + t148 * t166 + t151 * pkin(1) + 0; t144 * t154 + t146 * t156, t144 * t156 - t146 * t154, -(-t148 * t155 + t151 * t161) * t149 - t151 * t165, (-t153 * t159 + t170) * t157 + (t153 * t169 + t160) * t155 - t151 * t166 + t148 * pkin(1) + 0; t154 * t168 + (t154 * t163 + t155 * t156) * t150, t156 * t168 + (-t154 * t155 + t156 * t163) * t150, -t150 * t167 + t153 * t152, t147 * t153 + qJ(1) + 0 + (pkin(2) * t155 - pkin(10) * t167) * t150; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:22:30
	% EndTime: 2020-11-04 21:22:30
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (133->53), mult. (363->96), div. (0->0), fcn. (477->14), ass. (0->49)
	t188 = sin(pkin(7));
	t219 = pkin(10) * t188;
	t186 = sin(pkin(14));
	t218 = t186 * pkin(2);
	t190 = cos(pkin(14));
	t217 = t190 * pkin(2);
	t193 = cos(pkin(6));
	t216 = t188 * t193;
	t199 = cos(qJ(2));
	t215 = t188 * t199;
	t192 = cos(pkin(7));
	t185 = t192 * pkin(10) + pkin(9);
	t189 = sin(pkin(6));
	t214 = t189 * t185;
	t213 = t189 * t192;
	t196 = sin(qJ(2));
	t212 = t192 * t196;
	t211 = t192 * t199;
	t210 = t193 * t196;
	t209 = t193 * t199;
	t208 = t186 * t219;
	t207 = t190 * t219;
	t200 = -t188 * t189 + t192 * t209;
	t176 = -t200 * t186 - t190 * t212;
	t183 = t186 * t210 - t190 * t199;
	t195 = sin(qJ(3));
	t198 = cos(qJ(3));
	t173 = t176 * t198 + t195 * t183;
	t180 = (t186 * t209 + t190 * t196) * t188 + t186 * t213;
	t187 = sin(pkin(8));
	t191 = cos(pkin(8));
	t206 = t173 * t191 + t180 * t187;
	t205 = -t173 * t187 + t180 * t191;
	t177 = -t186 * t212 + t200 * t190;
	t184 = t186 * t199 + t190 * t210;
	t175 = t177 * t198 - t184 * t195;
	t181 = -(-t186 * t196 + t190 * t209) * t188 - t190 * t213;
	t204 = t175 * t191 + t181 * t187;
	t203 = -t175 * t187 + t181 * t191;
	t178 = t198 * t216 + (-t195 * t196 + t198 * t211) * t189;
	t182 = -t189 * t215 + t193 * t192;
	t202 = t178 * t191 + t182 * t187;
	t201 = -t178 * t187 + t182 * t191;
	t197 = cos(qJ(4));
	t194 = sin(qJ(4));
	t179 = t195 * t216 + (t195 * t211 + t196 * t198) * t189;
	t174 = t177 * t195 + t184 * t198;
	t172 = t176 * t195 - t198 * t183;
	t1 = [t172 * t197 + t206 * t194, -t172 * t194 + t206 * t197, t205, t172 * pkin(3) + (t193 * t208 + t217) * t199 + (-t193 * t218 + t207) * t196 + t186 * t214 + t190 * pkin(1) + 0 + t205 * pkin(11); t174 * t197 + t204 * t194, -t174 * t194 + t204 * t197, t203, t174 * pkin(3) + (-t193 * t207 + t218) * t199 + (t193 * t217 + t208) * t196 - t190 * t214 + t186 * pkin(1) + 0 + t203 * pkin(11); t179 * t197 + t202 * t194, -t179 * t194 + t202 * t197, t201, t179 * pkin(3) + t185 * t193 + qJ(1) + 0 + (pkin(2) * t196 - pkin(10) * t215) * t189 + t201 * pkin(11); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:22:30
	% EndTime: 2020-11-04 21:22:30
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (256->65), mult. (709->114), div. (0->0), fcn. (940->16), ass. (0->57)
	t245 = sin(pkin(7));
	t275 = pkin(10) * t245;
	t243 = sin(pkin(14));
	t274 = t243 * pkin(2);
	t247 = cos(pkin(14));
	t273 = t247 * pkin(2);
	t250 = cos(pkin(6));
	t272 = t245 * t250;
	t258 = cos(qJ(2));
	t271 = t245 * t258;
	t249 = cos(pkin(7));
	t242 = t249 * pkin(10) + pkin(9);
	t246 = sin(pkin(6));
	t270 = t246 * t242;
	t269 = t246 * t249;
	t254 = sin(qJ(2));
	t268 = t249 * t254;
	t267 = t249 * t258;
	t266 = t250 * t254;
	t265 = t250 * t258;
	t264 = t243 * t275;
	t263 = t247 * t275;
	t259 = -t245 * t246 + t249 * t265;
	t233 = -t259 * t243 - t247 * t268;
	t240 = t243 * t266 - t247 * t258;
	t253 = sin(qJ(3));
	t257 = cos(qJ(3));
	t229 = t233 * t257 + t253 * t240;
	t237 = (t243 * t265 + t247 * t254) * t245 + t243 * t269;
	t244 = sin(pkin(8));
	t248 = cos(pkin(8));
	t262 = t229 * t248 + t237 * t244;
	t226 = -t229 * t244 + t237 * t248;
	t234 = -t243 * t268 + t259 * t247;
	t241 = t243 * t258 + t247 * t266;
	t231 = t234 * t257 - t241 * t253;
	t238 = -(-t243 * t254 + t247 * t265) * t245 - t247 * t269;
	t261 = t231 * t248 + t238 * t244;
	t227 = -t231 * t244 + t238 * t248;
	t235 = t257 * t272 + (-t253 * t254 + t257 * t267) * t246;
	t239 = -t246 * t271 + t250 * t249;
	t260 = t235 * t248 + t239 * t244;
	t232 = -t235 * t244 + t239 * t248;
	t256 = cos(qJ(4));
	t255 = cos(qJ(5));
	t252 = sin(qJ(4));
	t251 = sin(qJ(5));
	t236 = t253 * t272 + (t253 * t267 + t254 * t257) * t246;
	t230 = t234 * t253 + t241 * t257;
	t228 = t233 * t253 - t257 * t240;
	t225 = t236 * t256 + t260 * t252;
	t224 = t236 * t252 - t260 * t256;
	t223 = t230 * t256 + t261 * t252;
	t222 = t230 * t252 - t261 * t256;
	t221 = t228 * t256 + t262 * t252;
	t220 = t228 * t252 - t262 * t256;
	t1 = [t221 * t255 + t226 * t251, -t221 * t251 + t226 * t255, t220, t221 * pkin(4) + t220 * pkin(12) + t228 * pkin(3) + (t250 * t264 + t273) * t258 + (-t250 * t274 + t263) * t254 + t243 * t270 + t247 * pkin(1) + 0 + t226 * pkin(11); t223 * t255 + t227 * t251, -t223 * t251 + t227 * t255, t222, t223 * pkin(4) + t222 * pkin(12) + t230 * pkin(3) + (-t250 * t263 + t274) * t258 + (t250 * t273 + t264) * t254 - t247 * t270 + t243 * pkin(1) + 0 + t227 * pkin(11); t225 * t255 + t232 * t251, -t225 * t251 + t232 * t255, t224, t236 * pkin(3) + t225 * pkin(4) + t224 * pkin(12) + t242 * t250 + qJ(1) + 0 + (pkin(2) * t254 - pkin(10) * t271) * t246 + t232 * pkin(11); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:22:30
	% EndTime: 2020-11-04 21:22:30
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (455->77), mult. (1269->132), div. (0->0), fcn. (1693->18), ass. (0->65)
	t307 = sin(pkin(7));
	t339 = pkin(10) * t307;
	t305 = sin(pkin(14));
	t338 = t305 * pkin(2);
	t309 = cos(pkin(14));
	t337 = t309 * pkin(2);
	t312 = cos(pkin(6));
	t336 = t307 * t312;
	t322 = cos(qJ(2));
	t335 = t307 * t322;
	t311 = cos(pkin(7));
	t304 = t311 * pkin(10) + pkin(9);
	t308 = sin(pkin(6));
	t334 = t308 * t304;
	t333 = t308 * t311;
	t317 = sin(qJ(2));
	t332 = t311 * t317;
	t331 = t311 * t322;
	t330 = t312 * t317;
	t329 = t312 * t322;
	t328 = t305 * t339;
	t327 = t309 * t339;
	t323 = -t307 * t308 + t311 * t329;
	t295 = -t323 * t305 - t309 * t332;
	t302 = t305 * t330 - t309 * t322;
	t316 = sin(qJ(3));
	t321 = cos(qJ(3));
	t291 = t295 * t321 + t316 * t302;
	t299 = (t305 * t329 + t309 * t317) * t307 + t305 * t333;
	t306 = sin(pkin(8));
	t310 = cos(pkin(8));
	t326 = t291 * t310 + t299 * t306;
	t288 = -t291 * t306 + t299 * t310;
	t296 = -t305 * t332 + t323 * t309;
	t303 = t305 * t322 + t309 * t330;
	t293 = t296 * t321 - t303 * t316;
	t300 = -(-t305 * t317 + t309 * t329) * t307 - t309 * t333;
	t325 = t293 * t310 + t300 * t306;
	t289 = -t293 * t306 + t300 * t310;
	t297 = t321 * t336 + (-t316 * t317 + t321 * t331) * t308;
	t301 = -t308 * t335 + t312 * t311;
	t324 = t297 * t310 + t301 * t306;
	t294 = -t297 * t306 + t301 * t310;
	t320 = cos(qJ(4));
	t319 = cos(qJ(5));
	t318 = cos(qJ(6));
	t315 = sin(qJ(4));
	t314 = sin(qJ(5));
	t313 = sin(qJ(6));
	t298 = t316 * t336 + (t316 * t331 + t317 * t321) * t308;
	t292 = t296 * t316 + t303 * t321;
	t290 = t295 * t316 - t321 * t302;
	t287 = t298 * t320 + t324 * t315;
	t286 = t298 * t315 - t324 * t320;
	t285 = t287 * t319 + t294 * t314;
	t284 = t287 * t314 - t294 * t319;
	t283 = t292 * t320 + t325 * t315;
	t282 = t292 * t315 - t325 * t320;
	t281 = t290 * t320 + t326 * t315;
	t280 = t290 * t315 - t326 * t320;
	t279 = t283 * t319 + t289 * t314;
	t278 = t283 * t314 - t289 * t319;
	t277 = t281 * t319 + t288 * t314;
	t276 = t281 * t314 - t288 * t319;
	t1 = [t277 * t318 + t280 * t313, -t277 * t313 + t280 * t318, t276, t277 * pkin(5) + t276 * pkin(13) + t281 * pkin(4) + t280 * pkin(12) + t290 * pkin(3) + (t312 * t328 + t337) * t322 + (-t312 * t338 + t327) * t317 + t305 * t334 + t309 * pkin(1) + 0 + t288 * pkin(11); t279 * t318 + t282 * t313, -t279 * t313 + t282 * t318, t278, t279 * pkin(5) + t278 * pkin(13) + t283 * pkin(4) + t282 * pkin(12) + t292 * pkin(3) + (-t312 * t327 + t338) * t322 + (t312 * t337 + t328) * t317 - t309 * t334 + t305 * pkin(1) + 0 + t289 * pkin(11); t285 * t318 + t286 * t313, -t285 * t313 + t286 * t318, t284, t298 * pkin(3) + t287 * pkin(4) + t285 * pkin(5) + t286 * pkin(12) + t284 * pkin(13) + t304 * t312 + qJ(1) + 0 + (pkin(2) * t317 - pkin(10) * t335) * t308 + t294 * pkin(11); 0, 0, 0, 1;];
	Tc_mdh = t1;
end