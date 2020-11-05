% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRR8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:48
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:31
	% EndTime: 2020-11-04 22:48:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:31
	% EndTime: 2020-11-04 22:48:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t132 = cos(qJ(1));
	t131 = sin(qJ(1));
	t1 = [t132, -t131, 0, 0; t131, t132, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:31
	% EndTime: 2020-11-04 22:48:31
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
	% StartTime: 2020-11-04 22:48:31
	% EndTime: 2020-11-04 22:48:32
	% DurationCPUTime: 0.07s
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
	% StartTime: 2020-11-04 22:48:32
	% EndTime: 2020-11-04 22:48:32
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (91->44), mult. (209->80), div. (0->0), fcn. (264->12), ass. (0->38)
	t189 = sin(pkin(7));
	t191 = cos(pkin(7));
	t194 = sin(qJ(3));
	t198 = cos(qJ(3));
	t204 = pkin(3) * t194 - pkin(11) * t198;
	t179 = t189 * pkin(10) - t191 * t204;
	t184 = pkin(3) * t198 + pkin(11) * t194 + pkin(2);
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
	t177 = -t190 * t212 + t192 * t202;
	t180 = t189 * t210 + t190 * t191;
	t193 = sin(qJ(4));
	t197 = cos(qJ(4));
	t203 = t177 * t193 + t180 * t197;
	t201 = t191 * pkin(10) + t189 * t204 + pkin(9);
	t196 = sin(qJ(1));
	t183 = t191 * t209 - t205;
	t182 = t191 * t210 - t213;
	t181 = -t191 * t192 + t199 * t213;
	t178 = t183 * t193 + t197 * t211;
	t176 = -t190 * t202 - t192 * t212;
	t175 = t179 * t195 + t184 * t199 + pkin(1);
	t174 = t190 * t201 - t192 * t216;
	t1 = [(-t177 * t196 - t183 * t200) * t197 + t193 * (t180 * t196 + t189 * t206), t200 * t178 + t196 * t203, (t182 * t196 + t191 * t206) * t198 + (-t192 * t195 * t196 + t199 * t200) * t194, t174 * t196 + t175 * t200 + 0; (t177 * t197 - t180 * t193) * t200 + t196 * (-t183 * t197 + t193 * t211), t196 * t178 - t200 * t203, (-t182 * t198 + t192 * t209) * t200 + t196 * (t191 * t207 + t208), -t174 * t200 + t175 * t196 + 0; -t176 * t197 - t181 * t193, t176 * t193 - t181 * t197, -t192 * t189 * t198 + (-t191 * t205 + t209) * t190, t190 * t216 + t201 * t192 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:32
	% EndTime: 2020-11-04 22:48:32
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (137->46), mult. (235->77), div. (0->0), fcn. (290->14), ass. (0->44)
	t263 = sin(qJ(4)) * pkin(4) + pkin(10);
	t235 = sin(pkin(7));
	t236 = sin(pkin(6));
	t261 = t235 * t236;
	t240 = sin(qJ(3));
	t260 = t235 * t240;
	t238 = cos(pkin(6));
	t244 = cos(qJ(2));
	t259 = t238 * t244;
	t241 = sin(qJ(2));
	t258 = t240 * t241;
	t257 = t240 * t244;
	t242 = sin(qJ(1));
	t256 = t241 * t242;
	t243 = cos(qJ(3));
	t255 = t241 * t243;
	t245 = cos(qJ(1));
	t254 = t241 * t245;
	t253 = t243 * t244;
	t237 = cos(pkin(7));
	t230 = cos(qJ(4)) * pkin(4) + pkin(3);
	t246 = pkin(11) + pkin(12);
	t249 = t230 * t240 - t246 * t243;
	t219 = t235 * t263 - t249 * t237;
	t227 = t230 * t243 + t246 * t240 + pkin(2);
	t252 = t219 * t244 - t227 * t241;
	t248 = t237 * t257 + t255;
	t221 = -t236 * t260 + t248 * t238;
	t228 = t237 * t258 - t253;
	t251 = t221 * t245 - t242 * t228;
	t250 = t221 * t242 + t245 * t228;
	t247 = t249 * t235 + t263 * t237 + pkin(9);
	t234 = qJ(4) + qJ(5);
	t232 = cos(t234);
	t231 = sin(t234);
	t226 = t237 * t259 - t261;
	t225 = -t238 * t237 + t244 * t261;
	t224 = t235 * t259 + t236 * t237;
	t223 = t224 * t245 - t235 * t256;
	t222 = t224 * t242 + t235 * t254;
	t220 = -t248 * t236 - t238 * t260;
	t218 = t219 * t241 + t227 * t244 + pkin(1);
	t217 = t236 * t247 + t252 * t238;
	t1 = [t222 * t231 - t250 * t232, t232 * t222 + t250 * t231, (t226 * t242 + t237 * t254) * t243 + (-t238 * t256 + t245 * t244) * t240, t217 * t242 + t218 * t245 + 0; -t223 * t231 + t251 * t232, -t232 * t223 - t251 * t231, (-t226 * t243 + t238 * t258) * t245 + t242 * (t237 * t255 + t257), -t217 * t245 + t218 * t242 + 0; -t220 * t232 - t231 * t225, t220 * t231 - t232 * t225, -t238 * t235 * t243 + (-t237 * t253 + t258) * t236, -t252 * t236 + t247 * t238 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:48:32
	% EndTime: 2020-11-04 22:48:32
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (231->61), mult. (438->95), div. (0->0), fcn. (563->16), ass. (0->55)
	t321 = sin(qJ(4)) * pkin(4) + pkin(10);
	t291 = sin(pkin(7));
	t292 = sin(pkin(6));
	t319 = t291 * t292;
	t297 = sin(qJ(3));
	t318 = t291 * t297;
	t294 = cos(pkin(6));
	t302 = cos(qJ(2));
	t317 = t294 * t302;
	t298 = sin(qJ(2));
	t316 = t297 * t298;
	t315 = t297 * t302;
	t299 = sin(qJ(1));
	t314 = t298 * t299;
	t301 = cos(qJ(3));
	t313 = t298 * t301;
	t303 = cos(qJ(1));
	t312 = t298 * t303;
	t311 = t301 * t302;
	t293 = cos(pkin(7));
	t286 = cos(qJ(4)) * pkin(4) + pkin(3);
	t304 = pkin(11) + pkin(12);
	t307 = t286 * t297 - t304 * t301;
	t274 = t291 * t321 - t307 * t293;
	t283 = t286 * t301 + t304 * t297 + pkin(2);
	t310 = t274 * t302 - t283 * t298;
	t306 = t293 * t315 + t313;
	t276 = -t292 * t318 + t306 * t294;
	t284 = t293 * t316 - t311;
	t309 = t276 * t303 - t299 * t284;
	t308 = t276 * t299 + t303 * t284;
	t305 = t307 * t291 + t321 * t293 + pkin(9);
	t300 = cos(qJ(6));
	t295 = sin(qJ(6));
	t290 = qJ(4) + qJ(5);
	t288 = cos(t290);
	t287 = sin(t290);
	t282 = t293 * t317 - t319;
	t281 = -t294 * t293 + t302 * t319;
	t280 = t291 * t317 + t292 * t293;
	t279 = t280 * t303 - t291 * t314;
	t278 = t280 * t299 + t291 * t312;
	t277 = -t294 * t291 * t301 + (-t293 * t311 + t316) * t292;
	t275 = -t306 * t292 - t294 * t318;
	t273 = (-t282 * t301 + t294 * t316) * t303 + t299 * (t293 * t313 + t315);
	t272 = (t282 * t299 + t293 * t312) * t301 + (-t294 * t314 + t303 * t302) * t297;
	t271 = t274 * t298 + t283 * t302 + pkin(1);
	t270 = -t275 * t288 - t287 * t281;
	t269 = t275 * t287 - t288 * t281;
	t268 = -t288 * t279 - t309 * t287;
	t267 = t278 * t287 - t308 * t288;
	t266 = -t279 * t287 + t309 * t288;
	t265 = t288 * t278 + t308 * t287;
	t264 = t292 * t305 + t310 * t294;
	t1 = [t267 * t300 + t272 * t295, -t267 * t295 + t272 * t300, -t265, t267 * pkin(5) - t265 * pkin(13) + t264 * t299 + t271 * t303 + 0; t266 * t300 + t273 * t295, -t266 * t295 + t273 * t300, -t268, t266 * pkin(5) - t268 * pkin(13) - t264 * t303 + t271 * t299 + 0; t270 * t300 + t277 * t295, -t270 * t295 + t277 * t300, -t269, t270 * pkin(5) - t269 * pkin(13) - t310 * t292 + t305 * t294 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end