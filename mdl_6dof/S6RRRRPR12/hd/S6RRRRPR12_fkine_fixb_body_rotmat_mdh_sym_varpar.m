% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPR12 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:41
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:01
	% EndTime: 2020-11-04 22:41:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:01
	% EndTime: 2020-11-04 22:41:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t132 = cos(qJ(1));
	t131 = sin(qJ(1));
	t1 = [t132, -t131, 0, 0; t131, t132, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:01
	% EndTime: 2020-11-04 22:41:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t133 = sin(pkin(6));
	t136 = sin(qJ(1));
	t144 = t136 * t133;
	t135 = sin(qJ(2));
	t143 = t136 * t135;
	t137 = cos(qJ(2));
	t142 = t136 * t137;
	t138 = cos(qJ(1));
	t141 = t138 * t133;
	t140 = t138 * t135;
	t139 = t138 * t137;
	t134 = cos(pkin(6));
	t1 = [-t134 * t143 + t139, -t134 * t142 - t140, t144, t138 * pkin(1) + pkin(9) * t144 + 0; t134 * t140 + t142, t134 * t139 - t143, -t141, t136 * pkin(1) - pkin(9) * t141 + 0; t133 * t135, t133 * t137, t134, t134 * pkin(9) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:01
	% EndTime: 2020-11-04 22:41:02
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (44->32), mult. (104->58), div. (0->0), fcn. (138->10), ass. (0->29)
	t151 = sin(pkin(7));
	t154 = cos(pkin(6));
	t173 = t151 * t154;
	t159 = cos(qJ(2));
	t172 = t151 * t159;
	t152 = sin(pkin(6));
	t153 = cos(pkin(7));
	t171 = t152 * t153;
	t170 = t154 * t159;
	t155 = sin(qJ(3));
	t156 = sin(qJ(2));
	t169 = t155 * t156;
	t168 = t155 * t159;
	t158 = cos(qJ(3));
	t167 = t156 * t158;
	t157 = sin(qJ(1));
	t166 = t157 * t156;
	t165 = t158 * t159;
	t160 = cos(qJ(1));
	t164 = t160 * t156;
	t163 = t160 * t159;
	t162 = -pkin(2) * t156 + pkin(10) * t172;
	t146 = -t151 * t152 + t153 * t170;
	t161 = t146 * t155 + t154 * t167;
	t150 = t153 * pkin(10) + pkin(9);
	t148 = t151 * t156 * pkin(10) + pkin(2) * t159 + pkin(1);
	t147 = -t153 * t169 + t165;
	t145 = t152 * t150 + t162 * t154;
	t1 = [t160 * t147 - t161 * t157, (-t146 * t157 - t153 * t164) * t158 - (-t154 * t166 + t163) * t155, (t157 * t170 + t164) * t151 + t157 * t171, t145 * t157 + t148 * t160 + 0; t157 * t147 + t161 * t160, (t146 * t158 - t154 * t169) * t160 - t157 * (t153 * t167 + t168), -(t154 * t163 - t166) * t151 - t160 * t171, -t145 * t160 + t148 * t157 + 0; t155 * t173 + (t153 * t168 + t167) * t152, t158 * t173 + (t153 * t165 - t169) * t152, -t152 * t172 + t154 * t153, t150 * t154 - t162 * t152 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:02
	% EndTime: 2020-11-04 22:41:02
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (91->44), mult. (209->80), div. (0->0), fcn. (264->12), ass. (0->38)
	t189 = sin(pkin(7));
	t191 = cos(pkin(7));
	t194 = sin(qJ(3));
	t198 = cos(qJ(3));
	t204 = pkin(3) * t194 - pkin(11) * t198;
	t179 = t189 * pkin(10) - t204 * t191;
	t184 = t198 * pkin(3) + pkin(11) * t194 + pkin(2);
	t195 = sin(qJ(2));
	t199 = cos(qJ(2));
	t216 = -t179 * t199 + t184 * t195;
	t190 = sin(pkin(6));
	t213 = t189 * t190;
	t212 = t189 * t194;
	t211 = t189 * t195;
	t192 = cos(pkin(6));
	t210 = t192 * t199;
	t209 = t194 * t195;
	t208 = t194 * t199;
	t207 = t195 * t198;
	t200 = cos(qJ(1));
	t206 = t195 * t200;
	t205 = t198 * t199;
	t202 = t191 * t208 + t207;
	t177 = -t190 * t212 + t202 * t192;
	t180 = t189 * t210 + t190 * t191;
	t193 = sin(qJ(4));
	t197 = cos(qJ(4));
	t203 = t177 * t193 + t197 * t180;
	t201 = t191 * pkin(10) + t204 * t189 + pkin(9);
	t196 = sin(qJ(1));
	t183 = t191 * t209 - t205;
	t182 = t191 * t210 - t213;
	t181 = -t192 * t191 + t199 * t213;
	t178 = t183 * t193 + t197 * t211;
	t176 = -t202 * t190 - t192 * t212;
	t175 = t179 * t195 + t184 * t199 + pkin(1);
	t174 = t190 * t201 - t216 * t192;
	t1 = [(-t177 * t196 - t200 * t183) * t197 + t193 * (t180 * t196 + t189 * t206), t178 * t200 + t203 * t196, (t182 * t196 + t191 * t206) * t198 + (-t196 * t192 * t195 + t200 * t199) * t194, t174 * t196 + t175 * t200 + 0; (t177 * t197 - t193 * t180) * t200 + t196 * (-t183 * t197 + t193 * t211), t196 * t178 - t203 * t200, (-t182 * t198 + t192 * t209) * t200 + t196 * (t191 * t207 + t208), -t174 * t200 + t175 * t196 + 0; -t176 * t197 - t193 * t181, t176 * t193 - t197 * t181, -t192 * t189 * t198 + (-t191 * t205 + t209) * t190, t216 * t190 + t201 * t192 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:02
	% EndTime: 2020-11-04 22:41:02
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (137->46), mult. (235->77), div. (0->0), fcn. (290->14), ass. (0->44)
	t234 = sin(pkin(7));
	t236 = cos(pkin(7));
	t230 = cos(qJ(4)) * pkin(4) + pkin(3);
	t238 = qJ(5) + pkin(11);
	t240 = sin(qJ(3));
	t243 = cos(qJ(3));
	t247 = t230 * t240 - t238 * t243;
	t251 = sin(qJ(4)) * pkin(4) + pkin(10);
	t261 = t247 * t234 + t251 * t236 + pkin(9);
	t235 = sin(pkin(6));
	t260 = t234 * t235;
	t259 = t234 * t240;
	t237 = cos(pkin(6));
	t244 = cos(qJ(2));
	t258 = t237 * t244;
	t241 = sin(qJ(2));
	t257 = t240 * t241;
	t256 = t240 * t244;
	t242 = sin(qJ(1));
	t255 = t241 * t242;
	t254 = t241 * t243;
	t245 = cos(qJ(1));
	t253 = t241 * t245;
	t252 = t243 * t244;
	t219 = t234 * t251 - t247 * t236;
	t224 = t230 * t243 + t238 * t240 + pkin(2);
	t250 = t219 * t244 - t224 * t241;
	t246 = t236 * t256 + t254;
	t221 = -t235 * t259 + t246 * t237;
	t228 = t236 * t257 - t252;
	t249 = t221 * t245 - t242 * t228;
	t248 = t221 * t242 + t245 * t228;
	t233 = qJ(4) + pkin(13);
	t232 = cos(t233);
	t231 = sin(t233);
	t227 = t236 * t258 - t260;
	t226 = -t237 * t236 + t244 * t260;
	t225 = t234 * t258 + t235 * t236;
	t223 = t225 * t245 - t234 * t255;
	t222 = t225 * t242 + t234 * t253;
	t220 = -t246 * t235 - t237 * t259;
	t218 = t219 * t241 + t224 * t244 + pkin(1);
	t217 = t261 * t235 + t250 * t237;
	t1 = [t231 * t222 - t248 * t232, t222 * t232 + t248 * t231, (t227 * t242 + t236 * t253) * t243 + (-t237 * t255 + t245 * t244) * t240, t217 * t242 + t218 * t245 + 0; -t223 * t231 + t249 * t232, -t223 * t232 - t249 * t231, (-t227 * t243 + t237 * t257) * t245 + t242 * (t236 * t254 + t256), -t217 * t245 + t218 * t242 + 0; -t220 * t232 - t231 * t226, t220 * t231 - t232 * t226, -t237 * t234 * t243 + (-t236 * t252 + t257) * t235, -t250 * t235 + t261 * t237 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:41:02
	% EndTime: 2020-11-04 22:41:02
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (231->61), mult. (438->95), div. (0->0), fcn. (563->16), ass. (0->55)
	t288 = sin(pkin(7));
	t290 = cos(pkin(7));
	t284 = cos(qJ(4)) * pkin(4) + pkin(3);
	t292 = qJ(5) + pkin(11);
	t295 = sin(qJ(3));
	t299 = cos(qJ(3));
	t303 = t284 * t295 - t292 * t299;
	t307 = sin(qJ(4)) * pkin(4) + pkin(10);
	t317 = t303 * t288 + t307 * t290 + pkin(9);
	t289 = sin(pkin(6));
	t316 = t288 * t289;
	t315 = t288 * t295;
	t291 = cos(pkin(6));
	t300 = cos(qJ(2));
	t314 = t291 * t300;
	t296 = sin(qJ(2));
	t313 = t295 * t296;
	t312 = t295 * t300;
	t297 = sin(qJ(1));
	t311 = t296 * t297;
	t310 = t296 * t299;
	t301 = cos(qJ(1));
	t309 = t296 * t301;
	t308 = t299 * t300;
	t272 = t288 * t307 - t303 * t290;
	t278 = t284 * t299 + t292 * t295 + pkin(2);
	t306 = t272 * t300 - t278 * t296;
	t302 = t290 * t312 + t310;
	t274 = -t289 * t315 + t302 * t291;
	t282 = t290 * t313 - t308;
	t305 = t274 * t301 - t297 * t282;
	t304 = t274 * t297 + t301 * t282;
	t298 = cos(qJ(6));
	t293 = sin(qJ(6));
	t287 = qJ(4) + pkin(13);
	t286 = cos(t287);
	t285 = sin(t287);
	t281 = t290 * t314 - t316;
	t280 = -t291 * t290 + t300 * t316;
	t279 = t288 * t314 + t289 * t290;
	t277 = t279 * t301 - t288 * t311;
	t276 = t279 * t297 + t288 * t309;
	t275 = -t291 * t288 * t299 + (-t290 * t308 + t313) * t289;
	t273 = -t302 * t289 - t291 * t315;
	t271 = (-t281 * t299 + t291 * t313) * t301 + t297 * (t290 * t310 + t312);
	t270 = (t281 * t297 + t290 * t309) * t299 + (-t291 * t311 + t301 * t300) * t295;
	t269 = t272 * t296 + t278 * t300 + pkin(1);
	t268 = -t273 * t286 - t285 * t280;
	t267 = t273 * t285 - t286 * t280;
	t266 = -t277 * t286 - t305 * t285;
	t265 = t285 * t276 - t304 * t286;
	t264 = -t277 * t285 + t305 * t286;
	t263 = t276 * t286 + t304 * t285;
	t262 = t317 * t289 + t306 * t291;
	t1 = [t265 * t298 + t270 * t293, -t265 * t293 + t270 * t298, -t263, t265 * pkin(5) - t263 * pkin(12) + t262 * t297 + t269 * t301 + 0; t264 * t298 + t271 * t293, -t264 * t293 + t271 * t298, -t266, t264 * pkin(5) - t266 * pkin(12) - t262 * t301 + t269 * t297 + 0; t268 * t298 + t275 * t293, -t268 * t293 + t275 * t298, -t267, t268 * pkin(5) - t267 * pkin(12) - t306 * t289 + t317 * t291 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end