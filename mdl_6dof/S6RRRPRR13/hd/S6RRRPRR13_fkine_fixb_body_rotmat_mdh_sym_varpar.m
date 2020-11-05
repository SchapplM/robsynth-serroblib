% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR13 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:33
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRR13_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR13_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:46
	% EndTime: 2020-11-04 22:33:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:46
	% EndTime: 2020-11-04 22:33:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t136 = cos(qJ(1));
	t135 = sin(qJ(1));
	t1 = [t136, -t135, 0, 0; t135, t136, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:46
	% EndTime: 2020-11-04 22:33:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t137 = sin(pkin(6));
	t140 = sin(qJ(1));
	t148 = t140 * t137;
	t139 = sin(qJ(2));
	t147 = t140 * t139;
	t141 = cos(qJ(2));
	t146 = t140 * t141;
	t142 = cos(qJ(1));
	t145 = t142 * t137;
	t144 = t142 * t139;
	t143 = t142 * t141;
	t138 = cos(pkin(6));
	t1 = [-t138 * t147 + t143, -t138 * t146 - t144, t148, t142 * pkin(1) + pkin(9) * t148 + 0; t138 * t144 + t146, t138 * t143 - t147, -t145, t140 * pkin(1) - pkin(9) * t145 + 0; t137 * t139, t137 * t141, t138, t138 * pkin(9) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:46
	% EndTime: 2020-11-04 22:33:46
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->32), mult. (104->58), div. (0->0), fcn. (138->10), ass. (0->29)
	t155 = sin(pkin(7));
	t158 = cos(pkin(6));
	t177 = t155 * t158;
	t163 = cos(qJ(2));
	t176 = t155 * t163;
	t156 = sin(pkin(6));
	t157 = cos(pkin(7));
	t175 = t156 * t157;
	t174 = t158 * t163;
	t159 = sin(qJ(3));
	t160 = sin(qJ(2));
	t173 = t159 * t160;
	t172 = t159 * t163;
	t162 = cos(qJ(3));
	t171 = t160 * t162;
	t161 = sin(qJ(1));
	t170 = t161 * t160;
	t169 = t162 * t163;
	t164 = cos(qJ(1));
	t168 = t164 * t160;
	t167 = t164 * t163;
	t166 = -pkin(2) * t160 + pkin(10) * t176;
	t150 = -t155 * t156 + t157 * t174;
	t165 = t150 * t159 + t158 * t171;
	t154 = t157 * pkin(10) + pkin(9);
	t152 = t155 * t160 * pkin(10) + pkin(2) * t163 + pkin(1);
	t151 = -t157 * t173 + t169;
	t149 = t156 * t154 + t166 * t158;
	t1 = [t164 * t151 - t165 * t161, (-t150 * t161 - t157 * t168) * t162 - (-t158 * t170 + t167) * t159, (t161 * t174 + t168) * t155 + t161 * t175, t149 * t161 + t152 * t164 + 0; t161 * t151 + t165 * t164, (t150 * t162 - t158 * t173) * t164 - t161 * (t157 * t171 + t172), -(t158 * t167 - t170) * t155 - t164 * t175, -t149 * t164 + t152 * t161 + 0; t159 * t177 + (t157 * t172 + t171) * t156, t162 * t177 + (t157 * t169 - t173) * t156, -t156 * t176 + t158 * t157, t154 * t158 - t166 * t156 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:46
	% EndTime: 2020-11-04 22:33:46
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (91->44), mult. (211->81), div. (0->0), fcn. (266->12), ass. (0->37)
	t191 = sin(pkin(7));
	t194 = cos(pkin(7));
	t196 = sin(qJ(3));
	t199 = cos(qJ(3));
	t204 = pkin(3) * t196 - qJ(4) * t199;
	t184 = t191 * pkin(10) - t204 * t194;
	t188 = pkin(3) * t199 + qJ(4) * t196 + pkin(2);
	t197 = sin(qJ(2));
	t200 = cos(qJ(2));
	t217 = -t184 * t200 + t188 * t197;
	t216 = t194 * pkin(10) + t204 * t191 + pkin(9);
	t190 = sin(pkin(13));
	t209 = t194 * t196;
	t193 = cos(pkin(13));
	t211 = t193 * t191;
	t185 = t190 * t209 + t211;
	t214 = t185 * t200;
	t212 = t191 * t196;
	t195 = cos(pkin(6));
	t210 = t194 * t195;
	t208 = t196 * t197;
	t207 = t197 * t199;
	t206 = t199 * t200;
	t187 = -t190 * t191 + t193 * t209;
	t203 = t187 * t200 + t193 * t207;
	t202 = t190 * t194 + t196 * t211;
	t201 = cos(qJ(1));
	t198 = sin(qJ(1));
	t192 = sin(pkin(6));
	t186 = -t191 * t192 + t200 * t210;
	t183 = -t197 * t187 + t193 * t206;
	t182 = t197 * t185 - t190 * t206;
	t181 = t184 * t197 + t188 * t200 + pkin(1);
	t180 = -t192 * t202 + t203 * t195;
	t179 = t192 * (-t190 * t212 + t193 * t194) + (t190 * t207 + t214) * t195;
	t178 = t192 * t216 - t217 * t195;
	t1 = [-t180 * t198 + t183 * t201, t179 * t198 + t182 * t201, (t201 * t197 * t194 + t186 * t198) * t199 + (-t198 * t195 * t197 + t201 * t200) * t196, t178 * t198 + t181 * t201 + 0; t180 * t201 + t198 * t183, -t179 * t201 + t198 * t182, (-t186 * t199 + t195 * t208) * t201 + t198 * (t194 * t207 + t196 * t200), -t178 * t201 + t181 * t198 + 0; t203 * t192 + t202 * t195, -t192 * t214 + t193 * t210 + (-t192 * t207 - t195 * t212) * t190, -t195 * t191 * t199 + (-t194 * t206 + t208) * t192, t217 * t192 + t216 * t195 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:46
	% EndTime: 2020-11-04 22:33:46
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (137->46), mult. (231->77), div. (0->0), fcn. (286->14), ass. (0->44)
	t231 = sin(pkin(13)) * pkin(4) + pkin(10);
	t236 = sin(pkin(7));
	t238 = cos(pkin(7));
	t232 = cos(pkin(13)) * pkin(4) + pkin(3);
	t240 = qJ(4) + pkin(11);
	t241 = sin(qJ(3));
	t244 = cos(qJ(3));
	t248 = t232 * t241 - t240 * t244;
	t262 = t231 * t238 + t248 * t236 + pkin(9);
	t237 = sin(pkin(6));
	t261 = t236 * t237;
	t260 = t236 * t241;
	t239 = cos(pkin(6));
	t245 = cos(qJ(2));
	t259 = t239 * t245;
	t242 = sin(qJ(2));
	t258 = t241 * t242;
	t257 = t241 * t245;
	t243 = sin(qJ(1));
	t256 = t242 * t243;
	t255 = t242 * t244;
	t246 = cos(qJ(1));
	t254 = t242 * t246;
	t253 = t244 * t245;
	t220 = t236 * t231 - t248 * t238;
	t225 = t232 * t244 + t240 * t241 + pkin(2);
	t251 = t220 * t245 - t225 * t242;
	t247 = t238 * t257 + t255;
	t222 = -t237 * t260 + t247 * t239;
	t229 = t238 * t258 - t253;
	t250 = t222 * t246 - t243 * t229;
	t249 = t222 * t243 + t246 * t229;
	t235 = pkin(13) + qJ(5);
	t234 = cos(t235);
	t233 = sin(t235);
	t228 = t238 * t259 - t261;
	t227 = -t239 * t238 + t245 * t261;
	t226 = t236 * t259 + t238 * t237;
	t224 = t226 * t246 - t236 * t256;
	t223 = t226 * t243 + t236 * t254;
	t221 = -t247 * t237 - t239 * t260;
	t219 = t220 * t242 + t225 * t245 + pkin(1);
	t218 = t262 * t237 + t251 * t239;
	t1 = [t223 * t233 - t249 * t234, t223 * t234 + t249 * t233, (t228 * t243 + t238 * t254) * t244 + (-t239 * t256 + t246 * t245) * t241, t218 * t243 + t219 * t246 + 0; -t224 * t233 + t250 * t234, -t224 * t234 - t250 * t233, (-t228 * t244 + t239 * t258) * t246 + t243 * (t238 * t255 + t257), -t218 * t246 + t219 * t243 + 0; -t221 * t234 - t233 * t227, t221 * t233 - t234 * t227, -t239 * t236 * t244 + (-t238 * t253 + t258) * t237, -t251 * t237 + t262 * t239 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:33:46
	% EndTime: 2020-11-04 22:33:46
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (231->61), mult. (434->95), div. (0->0), fcn. (559->16), ass. (0->55)
	t285 = sin(pkin(13)) * pkin(4) + pkin(10);
	t290 = sin(pkin(7));
	t292 = cos(pkin(7));
	t286 = cos(pkin(13)) * pkin(4) + pkin(3);
	t294 = qJ(4) + pkin(11);
	t296 = sin(qJ(3));
	t300 = cos(qJ(3));
	t304 = t286 * t296 - t294 * t300;
	t318 = t285 * t292 + t304 * t290 + pkin(9);
	t291 = sin(pkin(6));
	t317 = t290 * t291;
	t316 = t290 * t296;
	t293 = cos(pkin(6));
	t301 = cos(qJ(2));
	t315 = t293 * t301;
	t297 = sin(qJ(2));
	t314 = t296 * t297;
	t313 = t296 * t301;
	t298 = sin(qJ(1));
	t312 = t297 * t298;
	t311 = t297 * t300;
	t302 = cos(qJ(1));
	t310 = t297 * t302;
	t309 = t300 * t301;
	t273 = t290 * t285 - t304 * t292;
	t279 = t286 * t300 + t294 * t296 + pkin(2);
	t307 = t273 * t301 - t279 * t297;
	t303 = t292 * t313 + t311;
	t275 = -t291 * t316 + t303 * t293;
	t283 = t292 * t314 - t309;
	t306 = t275 * t302 - t298 * t283;
	t305 = t275 * t298 + t302 * t283;
	t299 = cos(qJ(6));
	t295 = sin(qJ(6));
	t289 = pkin(13) + qJ(5);
	t288 = cos(t289);
	t287 = sin(t289);
	t282 = t292 * t315 - t317;
	t281 = -t293 * t292 + t301 * t317;
	t280 = t290 * t315 + t292 * t291;
	t278 = t280 * t302 - t290 * t312;
	t277 = t280 * t298 + t290 * t310;
	t276 = -t293 * t290 * t300 + (-t292 * t309 + t314) * t291;
	t274 = -t303 * t291 - t293 * t316;
	t272 = (-t282 * t300 + t293 * t314) * t302 + t298 * (t292 * t311 + t313);
	t271 = (t282 * t298 + t292 * t310) * t300 + (-t293 * t312 + t302 * t301) * t296;
	t270 = t273 * t297 + t279 * t301 + pkin(1);
	t269 = -t274 * t288 - t287 * t281;
	t268 = t274 * t287 - t288 * t281;
	t267 = -t278 * t288 - t306 * t287;
	t266 = t277 * t287 - t305 * t288;
	t265 = -t278 * t287 + t306 * t288;
	t264 = t277 * t288 + t305 * t287;
	t263 = t318 * t291 + t307 * t293;
	t1 = [t266 * t299 + t271 * t295, -t266 * t295 + t271 * t299, -t264, t266 * pkin(5) - t264 * pkin(12) + t263 * t298 + t270 * t302 + 0; t265 * t299 + t272 * t295, -t265 * t295 + t272 * t299, -t267, t265 * pkin(5) - t267 * pkin(12) - t263 * t302 + t270 * t298 + 0; t269 * t299 + t276 * t295, -t269 * t295 + t276 * t299, -t268, t269 * pkin(5) - t268 * pkin(12) - t307 * t291 + t318 * t293 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end