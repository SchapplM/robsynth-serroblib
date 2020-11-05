% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPRRRR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:56
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PPRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:56:12
	% EndTime: 2020-11-04 20:56:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:56:12
	% EndTime: 2020-11-04 20:56:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t120 = cos(pkin(12));
	t119 = sin(pkin(12));
	t1 = [t120, -t119, 0, 0; t119, t120, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:56:12
	% EndTime: 2020-11-04 20:56:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t122 = sin(pkin(12));
	t123 = sin(pkin(6));
	t130 = t122 * t123;
	t126 = cos(pkin(6));
	t129 = t122 * t126;
	t125 = cos(pkin(12));
	t128 = t125 * t123;
	t127 = t125 * t126;
	t124 = cos(pkin(13));
	t121 = sin(pkin(13));
	t1 = [-t121 * t129 + t125 * t124, -t125 * t121 - t124 * t129, t130, t125 * pkin(1) + qJ(2) * t130 + 0; t121 * t127 + t122 * t124, -t122 * t121 + t124 * t127, -t128, t122 * pkin(1) - qJ(2) * t128 + 0; t123 * t121, t123 * t124, t126, t126 * qJ(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:56:12
	% EndTime: 2020-11-04 20:56:12
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (44->30), mult. (101->54), div. (0->0), fcn. (135->10), ass. (0->29)
	t142 = sin(pkin(7));
	t159 = pkin(8) * t142;
	t140 = sin(pkin(13));
	t145 = cos(pkin(12));
	t158 = t140 * t145;
	t141 = sin(pkin(12));
	t157 = t141 * t140;
	t147 = cos(pkin(6));
	t156 = t142 * t147;
	t143 = sin(pkin(6));
	t155 = t143 * t142;
	t144 = cos(pkin(13));
	t146 = cos(pkin(7));
	t154 = t144 * t146;
	t153 = t145 * t147;
	t152 = t146 * t143;
	t151 = t147 * t146;
	t137 = -t140 * pkin(2) + t144 * t159;
	t138 = t146 * pkin(8) + qJ(2);
	t150 = t137 * t147 + t143 * t138;
	t149 = cos(qJ(3));
	t148 = sin(qJ(3));
	t136 = t144 * pkin(2) + t140 * t159 + pkin(1);
	t135 = t145 * t144 - t147 * t157;
	t134 = t140 * t153 + t141 * t144;
	t133 = t144 * t151 - t155;
	t132 = -t133 * t141 - t146 * t158;
	t131 = t133 * t145 - t146 * t157;
	t1 = [t132 * t148 + t135 * t149, t132 * t149 - t135 * t148, (t144 * t156 + t152) * t141 + t142 * t158, t136 * t145 + t150 * t141 + 0; t131 * t148 + t134 * t149, t131 * t149 - t134 * t148, (-t144 * t153 + t157) * t142 - t145 * t152, t136 * t141 - t150 * t145 + 0; t148 * t156 + (t140 * t149 + t148 * t154) * t143, t149 * t156 + (-t140 * t148 + t149 * t154) * t143, -t144 * t155 + t151, -t137 * t143 + t138 * t147 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:56:12
	% EndTime: 2020-11-04 20:56:12
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (91->53), mult. (217->95), div. (0->0), fcn. (272->12), ass. (0->43)
	t180 = sin(pkin(13));
	t182 = sin(pkin(7));
	t205 = t180 * t182;
	t186 = cos(pkin(7));
	t204 = t180 * t186;
	t187 = cos(pkin(6));
	t203 = t180 * t187;
	t183 = sin(pkin(6));
	t202 = t183 * t182;
	t184 = cos(pkin(13));
	t201 = t184 * t186;
	t200 = t186 * t183;
	t199 = t187 * t182;
	t198 = t187 * t186;
	t166 = t184 * t198 - t202;
	t181 = sin(pkin(12));
	t185 = cos(pkin(12));
	t162 = t166 * t181 + t185 * t204;
	t169 = -t181 * t203 + t184 * t185;
	t189 = sin(qJ(3));
	t191 = cos(qJ(3));
	t197 = t162 * t189 - t169 * t191;
	t163 = -t166 * t185 + t181 * t204;
	t168 = t181 * t184 + t185 * t203;
	t196 = t163 * t189 - t168 * t191;
	t171 = pkin(8) * t182 * t184 - t180 * pkin(2);
	t177 = pkin(8) * t186 + qJ(2);
	t195 = t171 * t187 + t183 * t177;
	t175 = pkin(3) * t201 + pkin(9) * t180;
	t194 = pkin(3) * t202 - t175 * t187;
	t174 = -t180 * pkin(3) + pkin(9) * t201;
	t193 = pkin(9) * t202 - t174 * t187;
	t192 = t183 * t180 * t191 + (t184 * t200 + t199) * t189;
	t190 = cos(qJ(4));
	t188 = sin(qJ(4));
	t173 = pkin(3) * t204 - pkin(9) * t184;
	t172 = pkin(3) * t184 + pkin(9) * t204;
	t170 = pkin(2) * t184 + pkin(8) * t205 + pkin(1);
	t165 = t184 * t202 - t198;
	t164 = t184 * t199 + t200;
	t161 = t164 * t185 - t181 * t205;
	t160 = t164 * t181 + t185 * t205;
	t1 = [t160 * t188 - t190 * t197, t160 * t190 + t188 * t197, t162 * t191 + t169 * t189, (t185 * t172 - t181 * t193) * t191 + (-t185 * t173 + t181 * t194) * t189 + t195 * t181 + t170 * t185 + 0; -t188 * t161 - t190 * t196, -t190 * t161 + t188 * t196, t163 * t191 + t168 * t189, (t181 * t172 + t185 * t193) * t191 + (-t181 * t173 - t185 * t194) * t189 - t195 * t185 + t170 * t181 + 0; -t188 * t165 + t190 * t192, -t190 * t165 - t188 * t192, -t191 * t199 + (t180 * t189 - t191 * t201) * t183, (-pkin(9) * t199 - t174 * t183) * t191 + (pkin(3) * t199 + t175 * t183) * t189 - t171 * t183 + t177 * t187 + 0 + qJ(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:56:12
	% EndTime: 2020-11-04 20:56:12
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (167->68), mult. (411->113), div. (0->0), fcn. (536->14), ass. (0->54)
	t235 = sin(pkin(13));
	t237 = sin(pkin(7));
	t262 = t235 * t237;
	t241 = cos(pkin(7));
	t261 = t235 * t241;
	t242 = cos(pkin(6));
	t260 = t235 * t242;
	t238 = sin(pkin(6));
	t259 = t238 * t237;
	t239 = cos(pkin(13));
	t258 = t239 * t241;
	t257 = t241 * t238;
	t256 = t242 * t237;
	t255 = t242 * t241;
	t221 = t239 * t255 - t259;
	t236 = sin(pkin(12));
	t240 = cos(pkin(12));
	t217 = t221 * t236 + t240 * t261;
	t224 = -t236 * t260 + t240 * t239;
	t245 = sin(qJ(3));
	t248 = cos(qJ(3));
	t254 = t217 * t245 - t224 * t248;
	t218 = -t221 * t240 + t236 * t261;
	t223 = t236 * t239 + t240 * t260;
	t253 = t218 * t245 - t223 * t248;
	t226 = t239 * t237 * pkin(8) - t235 * pkin(2);
	t232 = t241 * pkin(8) + qJ(2);
	t252 = t226 * t242 + t238 * t232;
	t230 = pkin(3) * t258 + t235 * pkin(9);
	t251 = pkin(3) * t259 - t230 * t242;
	t229 = -t235 * pkin(3) + pkin(9) * t258;
	t250 = pkin(9) * t259 - t229 * t242;
	t249 = t238 * t235 * t248 + (t239 * t257 + t256) * t245;
	t247 = cos(qJ(4));
	t246 = cos(qJ(5));
	t244 = sin(qJ(4));
	t243 = sin(qJ(5));
	t228 = pkin(3) * t261 - t239 * pkin(9);
	t227 = t239 * pkin(3) + pkin(9) * t261;
	t225 = t239 * pkin(2) + pkin(8) * t262 + pkin(1);
	t220 = t239 * t259 - t255;
	t219 = t239 * t256 + t257;
	t216 = t219 * t240 - t236 * t262;
	t215 = t219 * t236 + t240 * t262;
	t214 = -t248 * t256 + (t235 * t245 - t248 * t258) * t238;
	t213 = -t247 * t220 - t249 * t244;
	t212 = -t244 * t220 + t249 * t247;
	t211 = t218 * t248 + t223 * t245;
	t210 = t217 * t248 + t224 * t245;
	t209 = -t247 * t216 + t253 * t244;
	t208 = t215 * t244 - t254 * t247;
	t207 = -t244 * t216 - t253 * t247;
	t206 = t215 * t247 + t254 * t244;
	t1 = [t208 * t246 + t210 * t243, -t208 * t243 + t210 * t246, -t206, t208 * pkin(4) - t206 * pkin(10) + (t240 * t227 - t250 * t236) * t248 + (-t240 * t228 + t251 * t236) * t245 + t252 * t236 + t225 * t240 + 0; t207 * t246 + t211 * t243, -t207 * t243 + t211 * t246, -t209, t207 * pkin(4) - t209 * pkin(10) + (t236 * t227 + t250 * t240) * t248 + (-t236 * t228 - t251 * t240) * t245 - t252 * t240 + t225 * t236 + 0; t212 * t246 + t214 * t243, -t212 * t243 + t214 * t246, -t213, t212 * pkin(4) - t213 * pkin(10) + (-pkin(9) * t256 - t229 * t238) * t248 + (pkin(3) * t256 + t230 * t238) * t245 - t226 * t238 + t232 * t242 + 0 + qJ(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:56:12
	% EndTime: 2020-11-04 20:56:12
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (198->74), mult. (449->118), div. (0->0), fcn. (584->16), ass. (0->58)
	t324 = pkin(5) * sin(qJ(5));
	t296 = sin(pkin(13));
	t298 = sin(pkin(7));
	t323 = t296 * t298;
	t302 = cos(pkin(7));
	t322 = t296 * t302;
	t303 = cos(pkin(6));
	t321 = t296 * t303;
	t299 = sin(pkin(6));
	t320 = t299 * t298;
	t300 = cos(pkin(13));
	t319 = t300 * t302;
	t318 = t302 * t299;
	t317 = t303 * t298;
	t316 = t303 * t302;
	t278 = t300 * t316 - t320;
	t297 = sin(pkin(12));
	t301 = cos(pkin(12));
	t274 = t278 * t297 + t301 * t322;
	t281 = -t297 * t321 + t301 * t300;
	t306 = sin(qJ(3));
	t308 = cos(qJ(3));
	t315 = t274 * t306 - t281 * t308;
	t275 = -t278 * t301 + t297 * t322;
	t280 = t297 * t300 + t301 * t321;
	t314 = t275 * t306 - t280 * t308;
	t283 = t300 * t298 * pkin(8) - t296 * pkin(2);
	t289 = t302 * pkin(8) + qJ(2);
	t313 = t283 * t303 + t299 * t289;
	t287 = pkin(3) * t319 + t296 * pkin(9);
	t312 = pkin(3) * t320 - t287 * t303;
	t286 = -t296 * pkin(3) + pkin(9) * t319;
	t311 = pkin(9) * t320 - t286 * t303;
	t310 = t299 * t296 * t308 + (t300 * t318 + t317) * t306;
	t309 = -pkin(11) - pkin(10);
	t307 = cos(qJ(4));
	t305 = sin(qJ(4));
	t295 = qJ(5) + qJ(6);
	t292 = cos(t295);
	t291 = sin(t295);
	t290 = cos(qJ(5)) * pkin(5) + pkin(4);
	t285 = pkin(3) * t322 - t300 * pkin(9);
	t284 = t300 * pkin(3) + pkin(9) * t322;
	t282 = t300 * pkin(2) + pkin(8) * t323 + pkin(1);
	t277 = t300 * t320 - t316;
	t276 = t300 * t317 + t318;
	t273 = t276 * t301 - t297 * t323;
	t272 = t276 * t297 + t301 * t323;
	t271 = -t308 * t317 + (t296 * t306 - t308 * t319) * t299;
	t270 = -t307 * t277 - t310 * t305;
	t269 = -t305 * t277 + t310 * t307;
	t268 = t275 * t308 + t280 * t306;
	t267 = t274 * t308 + t281 * t306;
	t266 = -t307 * t273 + t314 * t305;
	t265 = t272 * t305 - t315 * t307;
	t264 = -t305 * t273 - t314 * t307;
	t263 = t272 * t307 + t315 * t305;
	t1 = [t265 * t292 + t267 * t291, -t265 * t291 + t267 * t292, -t263, t265 * t290 + t263 * t309 + t267 * t324 + (t301 * t284 - t311 * t297) * t308 + (-t301 * t285 + t312 * t297) * t306 + t313 * t297 + t282 * t301 + 0; t264 * t292 + t268 * t291, -t264 * t291 + t268 * t292, -t266, t264 * t290 + t266 * t309 + t268 * t324 + (t297 * t284 + t311 * t301) * t308 + (-t297 * t285 - t312 * t301) * t306 - t313 * t301 + t282 * t297 + 0; t269 * t292 + t271 * t291, -t269 * t291 + t271 * t292, -t270, t269 * t290 + t270 * t309 + t271 * t324 + (-pkin(9) * t317 - t286 * t299) * t308 + (pkin(3) * t317 + t287 * t299) * t306 - t283 * t299 + t289 * t303 + 0 + qJ(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
end