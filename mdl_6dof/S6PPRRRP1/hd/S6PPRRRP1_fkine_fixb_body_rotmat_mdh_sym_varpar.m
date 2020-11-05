% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPRRRP1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:54
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PPRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:54:44
	% EndTime: 2020-11-04 20:54:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:54:44
	% EndTime: 2020-11-04 20:54:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t121 = cos(pkin(11));
	t120 = sin(pkin(11));
	t1 = [t121, -t120, 0, 0; t120, t121, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:54:44
	% EndTime: 2020-11-04 20:54:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t123 = sin(pkin(11));
	t124 = sin(pkin(6));
	t131 = t123 * t124;
	t127 = cos(pkin(6));
	t130 = t123 * t127;
	t126 = cos(pkin(11));
	t129 = t126 * t124;
	t128 = t126 * t127;
	t125 = cos(pkin(12));
	t122 = sin(pkin(12));
	t1 = [-t122 * t130 + t126 * t125, -t126 * t122 - t125 * t130, t131, t126 * pkin(1) + qJ(2) * t131 + 0; t122 * t128 + t123 * t125, -t123 * t122 + t125 * t128, -t129, t123 * pkin(1) - qJ(2) * t129 + 0; t124 * t122, t124 * t125, t127, t127 * qJ(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:54:44
	% EndTime: 2020-11-04 20:54:44
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (44->30), mult. (101->54), div. (0->0), fcn. (135->10), ass. (0->29)
	t143 = sin(pkin(7));
	t160 = pkin(8) * t143;
	t141 = sin(pkin(12));
	t146 = cos(pkin(11));
	t159 = t141 * t146;
	t142 = sin(pkin(11));
	t158 = t142 * t141;
	t148 = cos(pkin(6));
	t157 = t143 * t148;
	t144 = sin(pkin(6));
	t156 = t144 * t143;
	t145 = cos(pkin(12));
	t147 = cos(pkin(7));
	t155 = t145 * t147;
	t154 = t146 * t148;
	t153 = t147 * t144;
	t152 = t148 * t147;
	t138 = -t141 * pkin(2) + t145 * t160;
	t139 = t147 * pkin(8) + qJ(2);
	t151 = t138 * t148 + t144 * t139;
	t150 = cos(qJ(3));
	t149 = sin(qJ(3));
	t137 = t145 * pkin(2) + t141 * t160 + pkin(1);
	t136 = t146 * t145 - t148 * t158;
	t135 = t141 * t154 + t142 * t145;
	t134 = t145 * t152 - t156;
	t133 = -t134 * t142 - t147 * t159;
	t132 = t134 * t146 - t147 * t158;
	t1 = [t133 * t149 + t136 * t150, t133 * t150 - t136 * t149, (t145 * t157 + t153) * t142 + t143 * t159, t137 * t146 + t151 * t142 + 0; t132 * t149 + t135 * t150, t132 * t150 - t135 * t149, (-t145 * t154 + t158) * t143 - t146 * t153, t137 * t142 - t151 * t146 + 0; t149 * t157 + (t141 * t150 + t149 * t155) * t144, t150 * t157 + (-t141 * t149 + t150 * t155) * t144, -t145 * t156 + t152, -t138 * t144 + t139 * t148 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:54:44
	% EndTime: 2020-11-04 20:54:45
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (91->53), mult. (217->95), div. (0->0), fcn. (272->12), ass. (0->43)
	t181 = sin(pkin(12));
	t183 = sin(pkin(7));
	t206 = t181 * t183;
	t187 = cos(pkin(7));
	t205 = t181 * t187;
	t188 = cos(pkin(6));
	t204 = t181 * t188;
	t184 = sin(pkin(6));
	t203 = t184 * t183;
	t185 = cos(pkin(12));
	t202 = t185 * t187;
	t201 = t187 * t184;
	t200 = t188 * t183;
	t199 = t188 * t187;
	t168 = t185 * t199 - t203;
	t182 = sin(pkin(11));
	t186 = cos(pkin(11));
	t163 = t168 * t182 + t186 * t205;
	t170 = -t182 * t204 + t186 * t185;
	t190 = sin(qJ(3));
	t192 = cos(qJ(3));
	t198 = t163 * t190 - t170 * t192;
	t164 = -t168 * t186 + t182 * t205;
	t169 = t182 * t185 + t186 * t204;
	t197 = t164 * t190 - t169 * t192;
	t172 = t185 * t183 * pkin(8) - t181 * pkin(2);
	t178 = t187 * pkin(8) + qJ(2);
	t196 = t172 * t188 + t184 * t178;
	t176 = pkin(3) * t202 + t181 * pkin(9);
	t195 = pkin(3) * t203 - t176 * t188;
	t175 = -t181 * pkin(3) + pkin(9) * t202;
	t194 = pkin(9) * t203 - t175 * t188;
	t193 = t184 * t181 * t192 + (t185 * t201 + t200) * t190;
	t191 = cos(qJ(4));
	t189 = sin(qJ(4));
	t174 = pkin(3) * t205 - t185 * pkin(9);
	t173 = t185 * pkin(3) + pkin(9) * t205;
	t171 = t185 * pkin(2) + pkin(8) * t206 + pkin(1);
	t166 = t185 * t203 - t199;
	t165 = t185 * t200 + t201;
	t162 = t165 * t186 - t182 * t206;
	t161 = t165 * t182 + t186 * t206;
	t1 = [t161 * t189 - t198 * t191, t161 * t191 + t198 * t189, t163 * t192 + t170 * t190, (t186 * t173 - t194 * t182) * t192 + (-t186 * t174 + t195 * t182) * t190 + t196 * t182 + t171 * t186 + 0; -t189 * t162 - t197 * t191, -t162 * t191 + t197 * t189, t164 * t192 + t169 * t190, (t182 * t173 + t194 * t186) * t192 + (-t182 * t174 - t195 * t186) * t190 - t196 * t186 + t171 * t182 + 0; -t189 * t166 + t193 * t191, -t191 * t166 - t193 * t189, -t192 * t200 + (t181 * t190 - t192 * t202) * t184, (-pkin(9) * t200 - t175 * t184) * t192 + (pkin(3) * t200 + t176 * t184) * t190 - t172 * t184 + t178 * t188 + 0 + qJ(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:54:45
	% EndTime: 2020-11-04 20:54:45
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (167->68), mult. (411->113), div. (0->0), fcn. (536->14), ass. (0->53)
	t236 = sin(pkin(12));
	t238 = sin(pkin(7));
	t262 = t236 * t238;
	t242 = cos(pkin(7));
	t261 = t236 * t242;
	t243 = cos(pkin(6));
	t260 = t236 * t243;
	t259 = t238 * t243;
	t239 = sin(pkin(6));
	t258 = t239 * t238;
	t240 = cos(pkin(12));
	t257 = t240 * t242;
	t256 = t242 * t243;
	t223 = t240 * t256 - t258;
	t237 = sin(pkin(11));
	t241 = cos(pkin(11));
	t218 = t223 * t237 + t241 * t261;
	t225 = -t237 * t260 + t240 * t241;
	t246 = sin(qJ(3));
	t249 = cos(qJ(3));
	t255 = t218 * t246 - t225 * t249;
	t219 = -t223 * t241 + t237 * t261;
	t224 = t237 * t240 + t241 * t260;
	t254 = t219 * t246 - t224 * t249;
	t227 = pkin(8) * t238 * t240 - t236 * pkin(2);
	t233 = pkin(8) * t242 + qJ(2);
	t253 = t227 * t243 + t239 * t233;
	t231 = pkin(3) * t257 + pkin(9) * t236;
	t252 = pkin(3) * t258 - t231 * t243;
	t230 = -t236 * pkin(3) + pkin(9) * t257;
	t251 = pkin(9) * t258 - t230 * t243;
	t250 = t236 * t239 * t249 + (t239 * t257 + t259) * t246;
	t248 = cos(qJ(4));
	t247 = cos(qJ(5));
	t245 = sin(qJ(4));
	t244 = sin(qJ(5));
	t229 = pkin(3) * t261 - pkin(9) * t240;
	t228 = pkin(3) * t240 + pkin(9) * t261;
	t226 = pkin(2) * t240 + pkin(8) * t262 + pkin(1);
	t221 = t240 * t258 - t256;
	t220 = t239 * t242 + t240 * t259;
	t217 = t220 * t241 - t237 * t262;
	t216 = t220 * t237 + t241 * t262;
	t215 = -t249 * t259 + (t236 * t246 - t249 * t257) * t239;
	t214 = -t248 * t221 - t250 * t245;
	t213 = -t245 * t221 + t250 * t248;
	t212 = t219 * t249 + t224 * t246;
	t211 = t218 * t249 + t225 * t246;
	t210 = -t217 * t248 + t254 * t245;
	t209 = t216 * t245 - t255 * t248;
	t208 = -t245 * t217 - t254 * t248;
	t207 = t216 * t248 + t255 * t245;
	t1 = [t209 * t247 + t211 * t244, -t209 * t244 + t211 * t247, -t207, t209 * pkin(4) - t207 * pkin(10) + (t241 * t228 - t251 * t237) * t249 + (-t241 * t229 + t252 * t237) * t246 + t253 * t237 + t226 * t241 + 0; t208 * t247 + t212 * t244, -t208 * t244 + t212 * t247, -t210, t208 * pkin(4) - t210 * pkin(10) + (t237 * t228 + t251 * t241) * t249 + (-t237 * t229 - t252 * t241) * t246 - t253 * t241 + t226 * t237 + 0; t213 * t247 + t215 * t244, -t213 * t244 + t215 * t247, -t214, t213 * pkin(4) - t214 * pkin(10) + (-pkin(9) * t259 - t230 * t239) * t249 + (pkin(3) * t259 + t231 * t239) * t246 - t227 * t239 + t233 * t243 + 0 + qJ(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:54:45
	% EndTime: 2020-11-04 20:54:45
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (186->73), mult. (449->117), div. (0->0), fcn. (584->14), ass. (0->59)
	t297 = cos(pkin(12));
	t299 = cos(pkin(7));
	t300 = cos(pkin(6));
	t314 = t300 * t299;
	t295 = sin(pkin(7));
	t296 = sin(pkin(6));
	t318 = t296 * t295;
	t279 = t297 * t314 - t318;
	t294 = sin(pkin(11));
	t298 = cos(pkin(11));
	t293 = sin(pkin(12));
	t320 = t293 * t299;
	t274 = t279 * t294 + t298 * t320;
	t319 = t293 * t300;
	t281 = -t294 * t319 + t298 * t297;
	t304 = sin(qJ(3));
	t307 = cos(qJ(3));
	t267 = t274 * t307 + t281 * t304;
	t302 = sin(qJ(5));
	t324 = t267 * t302;
	t275 = -t279 * t298 + t294 * t320;
	t280 = t294 * t297 + t298 * t319;
	t268 = t275 * t307 + t280 * t304;
	t323 = t268 * t302;
	t315 = t300 * t295;
	t317 = t297 * t299;
	t271 = -t307 * t315 + (t293 * t304 - t307 * t317) * t296;
	t322 = t271 * t302;
	t321 = t293 * t295;
	t316 = t299 * t296;
	t313 = t274 * t304 - t281 * t307;
	t312 = t275 * t304 - t280 * t307;
	t283 = t297 * t295 * pkin(8) - t293 * pkin(2);
	t289 = t299 * pkin(8) + qJ(2);
	t311 = t283 * t300 + t296 * t289;
	t287 = pkin(3) * t317 + t293 * pkin(9);
	t310 = pkin(3) * t318 - t287 * t300;
	t286 = -t293 * pkin(3) + pkin(9) * t317;
	t309 = pkin(9) * t318 - t286 * t300;
	t308 = t296 * t293 * t307 + (t297 * t316 + t315) * t304;
	t306 = cos(qJ(4));
	t305 = cos(qJ(5));
	t303 = sin(qJ(4));
	t301 = -qJ(6) - pkin(10);
	t290 = t305 * pkin(5) + pkin(4);
	t285 = pkin(3) * t320 - t297 * pkin(9);
	t284 = t297 * pkin(3) + pkin(9) * t320;
	t282 = t297 * pkin(2) + pkin(8) * t321 + pkin(1);
	t277 = t297 * t318 - t314;
	t276 = t297 * t315 + t316;
	t273 = t276 * t298 - t294 * t321;
	t272 = t276 * t294 + t298 * t321;
	t270 = -t306 * t277 - t308 * t303;
	t269 = -t303 * t277 + t308 * t306;
	t266 = -t273 * t306 + t312 * t303;
	t265 = t272 * t303 - t313 * t306;
	t264 = -t303 * t273 - t312 * t306;
	t263 = t272 * t306 + t313 * t303;
	t1 = [t265 * t305 + t324, -t265 * t302 + t267 * t305, -t263, t265 * t290 + t263 * t301 + pkin(5) * t324 + (t298 * t284 - t309 * t294) * t307 + (-t298 * t285 + t310 * t294) * t304 + t311 * t294 + t282 * t298 + 0; t264 * t305 + t323, -t264 * t302 + t268 * t305, -t266, t264 * t290 + t266 * t301 + pkin(5) * t323 + (t294 * t284 + t309 * t298) * t307 + (-t294 * t285 - t310 * t298) * t304 - t311 * t298 + t282 * t294 + 0; t269 * t305 + t322, -t269 * t302 + t271 * t305, -t270, t269 * t290 + t270 * t301 + pkin(5) * t322 + (-pkin(9) * t315 - t286 * t296) * t307 + (pkin(3) * t315 + t287 * t296) * t304 - t283 * t296 + t289 * t300 + 0 + qJ(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
end