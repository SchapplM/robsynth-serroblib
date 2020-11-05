% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR14 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:28
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRR14_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR14_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:54
	% EndTime: 2020-11-04 20:28:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:54
	% EndTime: 2020-11-04 20:28:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t151 = cos(qJ(1));
	t150 = sin(qJ(1));
	t1 = [t151, -t150, 0, 0; t150, t151, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:54
	% EndTime: 2020-11-04 20:28:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t152 = sin(pkin(11));
	t156 = sin(qJ(1));
	t163 = t156 * t152;
	t153 = sin(pkin(5));
	t162 = t156 * t153;
	t154 = cos(pkin(11));
	t161 = t156 * t154;
	t157 = cos(qJ(1));
	t160 = t157 * t152;
	t159 = t157 * t153;
	t158 = t157 * t154;
	t155 = cos(pkin(5));
	t1 = [-t155 * t163 + t158, -t155 * t161 - t160, t162, t157 * pkin(1) + qJ(2) * t162 + 0; t155 * t160 + t161, t155 * t158 - t163, -t159, t156 * pkin(1) - qJ(2) * t159 + 0; t153 * t152, t153 * t154, t155, t155 * qJ(2) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:54
	% EndTime: 2020-11-04 20:28:54
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (44->32), mult. (101->58), div. (0->0), fcn. (135->10), ass. (0->31)
	t172 = sin(pkin(6));
	t195 = pkin(8) * t172;
	t171 = sin(pkin(11));
	t177 = sin(qJ(3));
	t194 = t171 * t177;
	t179 = cos(qJ(3));
	t193 = t171 * t179;
	t173 = sin(pkin(5));
	t192 = t173 * t172;
	t175 = cos(pkin(6));
	t191 = t173 * t175;
	t176 = cos(pkin(5));
	t190 = t176 * t175;
	t189 = t176 * t177;
	t188 = t176 * t179;
	t174 = cos(pkin(11));
	t187 = t177 * t174;
	t178 = sin(qJ(1));
	t186 = t178 * t171;
	t185 = t179 * t174;
	t180 = cos(qJ(1));
	t184 = t180 * t171;
	t183 = t180 * t174;
	t167 = -t171 * pkin(2) + t174 * t195;
	t169 = t175 * pkin(8) + qJ(2);
	t182 = t167 * t176 + t173 * t169;
	t164 = t174 * t190 - t192;
	t181 = t164 * t177 + t171 * t188;
	t166 = t174 * pkin(2) + t171 * t195 + pkin(1);
	t165 = t175 * t194 - t185;
	t1 = [-t180 * t165 - t181 * t178, (-t164 * t178 - t175 * t184) * t179 + t177 * (t176 * t186 - t183), (t178 * t176 * t174 + t184) * t172 + t178 * t191, t166 * t180 + t182 * t178 + 0; -t178 * t165 + t181 * t180, (t164 * t179 - t171 * t189) * t180 - t178 * (t175 * t193 + t187), -(t176 * t183 - t186) * t172 - t180 * t191, t166 * t178 - t182 * t180 + 0; t172 * t189 + (t175 * t187 + t193) * t173, t172 * t188 + (t175 * t185 - t194) * t173, -t174 * t192 + t190, -t167 * t173 + t169 * t176 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:54
	% EndTime: 2020-11-04 20:28:54
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (91->53), mult. (213->97), div. (0->0), fcn. (268->12), ass. (0->38)
	t215 = cos(pkin(11));
	t216 = cos(pkin(6));
	t217 = cos(pkin(5));
	t231 = t217 * t216;
	t213 = sin(pkin(6));
	t214 = sin(pkin(5));
	t234 = t214 * t213;
	t203 = t215 * t231 - t234;
	t212 = sin(pkin(11));
	t219 = sin(qJ(3));
	t222 = cos(qJ(3));
	t230 = t217 * t222;
	t197 = t203 * t219 + t212 * t230;
	t232 = t217 * t213;
	t202 = t214 * t216 + t215 * t232;
	t218 = sin(qJ(4));
	t221 = cos(qJ(4));
	t241 = t197 * t218 + t221 * t202;
	t205 = t215 * t213 * pkin(8) - t212 * pkin(2);
	t233 = t215 * t216;
	t206 = -t212 * pkin(3) + pkin(9) * t233;
	t207 = pkin(3) * t233 + t212 * pkin(9);
	t209 = t216 * pkin(8) + qJ(2);
	t240 = (pkin(3) * t234 - t207 * t217) * t219 - (pkin(9) * t234 - t206 * t217) * t222 + t205 * t217 + t214 * t209;
	t239 = t212 * t213;
	t238 = t212 * t216;
	t237 = t212 * t219;
	t236 = t212 * t222;
	t223 = cos(qJ(1));
	t235 = t212 * t223;
	t228 = t222 * t215;
	t224 = (t214 * t233 + t232) * t219 + t214 * t236;
	t220 = sin(qJ(1));
	t204 = -t216 * t237 + t228;
	t200 = t215 * t234 - t231;
	t199 = t204 * t218 - t221 * t239;
	t196 = (t215 * pkin(3) + pkin(9) * t238) * t222 + (-pkin(3) * t238 + t215 * pkin(9)) * t219 + pkin(8) * t239 + t215 * pkin(2) + pkin(1);
	t1 = [(-t197 * t220 + t223 * t204) * t221 + t218 * (t202 * t220 + t213 * t235), -t199 * t223 + t241 * t220, (t203 * t220 + t216 * t235) * t222 - t219 * (t220 * t217 * t212 - t223 * t215), t196 * t223 + t240 * t220 + 0; (t197 * t221 - t218 * t202) * t223 + t220 * (t204 * t221 + t218 * t239), -t220 * t199 - t241 * t223, (-t203 * t222 + t217 * t237) * t223 + t220 * (t219 * t215 + t216 * t236), t196 * t220 - t240 * t223 + 0; -t218 * t200 + t224 * t221, -t221 * t200 - t224 * t218, -t213 * t230 + (-t216 * t228 + t237) * t214, (-pkin(9) * t232 - t206 * t214) * t222 + (pkin(3) * t232 + t207 * t214) * t219 - t205 * t214 + t209 * t217 + 0 + pkin(7); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:28:54
	% EndTime: 2020-11-04 20:28:54
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (167->66), mult. (417->110), div. (0->0), fcn. (518->14), ass. (0->53)
	t267 = cos(pkin(11));
	t268 = cos(pkin(6));
	t269 = cos(pkin(5));
	t292 = t269 * t268;
	t265 = sin(pkin(6));
	t266 = sin(pkin(5));
	t297 = t266 * t265;
	t255 = t267 * t292 - t297;
	t261 = t268 * pkin(8) + qJ(2);
	t272 = sin(qJ(3));
	t264 = sin(pkin(11));
	t276 = cos(qJ(3));
	t295 = t267 * t268;
	t275 = cos(qJ(4));
	t299 = t264 * t275;
	t271 = sin(qJ(4));
	t302 = pkin(10) * t271;
	t279 = (pkin(4) * t299 - pkin(9) * t295 + (pkin(3) + t302) * t264) * t276 - t267 * t265 * pkin(8) + t264 * pkin(2) + (pkin(3) * t295 + t264 * pkin(9)) * t272;
	t287 = pkin(4) * t275 + t302;
	t293 = t269 * t265;
	t254 = t266 * t268 + t267 * t293;
	t290 = t275 * t254;
	t291 = t271 * t254;
	t296 = t266 * t276;
	t308 = -t265 * pkin(9) * t296 + pkin(4) * t291 - pkin(10) * t290 + t266 * t261 + (pkin(3) * t297 - t287 * t255) * t272 - t279 * t269;
	t298 = t264 * t276;
	t282 = t255 * t272 + t269 * t298;
	t306 = t282 * t271 + t290;
	t303 = pkin(9) * t276;
	t253 = t266 * t295 + t293;
	t301 = t253 * t272;
	t300 = t264 * t272;
	t294 = t268 * t272;
	t289 = t276 * t267;
	t286 = -pkin(4) * t271 + pkin(10) * t275;
	t243 = t282 * t275 - t291;
	t247 = t255 * t276 - t269 * t300;
	t270 = sin(qJ(5));
	t274 = cos(qJ(5));
	t285 = t243 * t270 + t247 * t274;
	t284 = pkin(3) + t287;
	t283 = t264 * t296 + t301;
	t277 = cos(qJ(1));
	t273 = sin(qJ(1));
	t256 = t272 * t267 + t268 * t298;
	t252 = t267 * t297 - t292;
	t249 = t275 * t289 - t264 * (-t265 * t271 + t275 * t294);
	t248 = (-t264 * t294 + t289) * t271 - t265 * t299;
	t246 = t253 * t276 - t266 * t300;
	t245 = -t249 * t270 + t274 * t256;
	t244 = -t271 * t252 + t283 * t275;
	t242 = pkin(1) + (pkin(9) * t272 + t284 * t276 + pkin(2)) * t267 + ((pkin(8) - t286) * t265 + (-t284 * t272 + t303) * t268) * t264;
	t1 = [(-t243 * t273 + t249 * t277) * t274 + t270 * (t247 * t273 + t277 * t256), t245 * t277 + t285 * t273, t248 * t277 - t306 * t273, t242 * t277 + t308 * t273 + 0; (t243 * t274 - t270 * t247) * t277 + (t249 * t274 + t270 * t256) * t273, t245 * t273 - t285 * t277, t273 * t248 + t306 * t277, t242 * t273 - t308 * t277 + 0; t244 * t274 - t246 * t270, -t244 * t270 - t246 * t274, t275 * t252 + t283 * t271, pkin(7) + 0 + t287 * t301 + (t261 + (pkin(3) * t272 - t303) * t265) * t269 + t286 * t252 + t279 * t266; 0, 0, 0, 1;];
	Tc_mdh = t1;
end