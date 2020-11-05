% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:12
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:12:13
	% EndTime: 2020-11-04 21:12:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:12:13
	% EndTime: 2020-11-04 21:12:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t111 = cos(pkin(12));
	t110 = sin(pkin(12));
	t1 = [t111, -t110, 0, 0; t110, t111, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:12:14
	% EndTime: 2020-11-04 21:12:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t112 = sin(pkin(12));
	t113 = sin(pkin(6));
	t121 = t112 * t113;
	t114 = cos(pkin(12));
	t120 = t114 * t113;
	t115 = cos(pkin(6));
	t116 = sin(qJ(2));
	t119 = t115 * t116;
	t117 = cos(qJ(2));
	t118 = t115 * t117;
	t1 = [-t112 * t119 + t114 * t117, -t112 * t118 - t114 * t116, t121, t114 * pkin(1) + pkin(8) * t121 + 0; t112 * t117 + t114 * t119, -t112 * t116 + t114 * t118, -t120, t112 * pkin(1) - pkin(8) * t120 + 0; t113 * t116, t113 * t117, t115, t115 * pkin(8) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:12:14
	% EndTime: 2020-11-04 21:12:14
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->35), mult. (112->66), div. (0->0), fcn. (146->10), ass. (0->30)
	t128 = sin(pkin(7));
	t150 = pkin(9) * t128;
	t127 = sin(pkin(12));
	t149 = t127 * pkin(2);
	t130 = cos(pkin(12));
	t148 = t130 * pkin(2);
	t132 = cos(pkin(6));
	t147 = t128 * t132;
	t136 = cos(qJ(2));
	t146 = t128 * t136;
	t131 = cos(pkin(7));
	t126 = t131 * pkin(9) + pkin(8);
	t129 = sin(pkin(6));
	t145 = t129 * t126;
	t144 = t129 * t131;
	t134 = sin(qJ(2));
	t143 = t131 * t134;
	t142 = t131 * t136;
	t141 = t132 * t134;
	t140 = t132 * t136;
	t139 = t127 * t150;
	t138 = t130 * t150;
	t137 = -t128 * t129 + t131 * t140;
	t135 = cos(qJ(3));
	t133 = sin(qJ(3));
	t125 = t127 * t136 + t130 * t141;
	t124 = t127 * t141 - t130 * t136;
	t123 = -t127 * t143 + t137 * t130;
	t122 = -t137 * t127 - t130 * t143;
	t1 = [t122 * t133 - t135 * t124, t122 * t135 + t133 * t124, (t127 * t140 + t130 * t134) * t128 + t127 * t144, (t132 * t139 + t148) * t136 + (-t132 * t149 + t138) * t134 + t127 * t145 + t130 * pkin(1) + 0; t123 * t133 + t125 * t135, t123 * t135 - t125 * t133, -(-t127 * t134 + t130 * t140) * t128 - t130 * t144, (-t132 * t138 + t149) * t136 + (t132 * t148 + t139) * t134 - t130 * t145 + t127 * pkin(1) + 0; t133 * t147 + (t133 * t142 + t134 * t135) * t129, t135 * t147 + (-t133 * t134 + t135 * t142) * t129, -t129 * t146 + t132 * t131, t126 * t132 + qJ(1) + 0 + (pkin(2) * t134 - pkin(9) * t146) * t129; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:12:14
	% EndTime: 2020-11-04 21:12:14
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (91->49), mult. (238->84), div. (0->0), fcn. (313->12), ass. (0->41)
	t167 = sin(pkin(7));
	t190 = pkin(9) * t167;
	t166 = sin(pkin(12));
	t189 = t166 * pkin(2);
	t170 = cos(pkin(12));
	t188 = t170 * pkin(2);
	t172 = cos(pkin(6));
	t187 = t167 * t172;
	t176 = cos(qJ(2));
	t186 = t167 * t176;
	t171 = cos(pkin(7));
	t164 = t171 * pkin(9) + pkin(8);
	t168 = sin(pkin(6));
	t185 = t168 * t164;
	t184 = t168 * t171;
	t174 = sin(qJ(2));
	t183 = t171 * t174;
	t182 = t171 * t176;
	t181 = t172 * t174;
	t180 = t172 * t176;
	t179 = t166 * t190;
	t178 = t170 * t190;
	t177 = -t167 * t168 + t171 * t180;
	t175 = cos(qJ(3));
	t173 = sin(qJ(3));
	t169 = cos(pkin(13));
	t165 = sin(pkin(13));
	t163 = t166 * t176 + t170 * t181;
	t162 = t166 * t181 - t170 * t176;
	t161 = -t168 * t186 + t172 * t171;
	t160 = -(-t166 * t174 + t170 * t180) * t167 - t170 * t184;
	t159 = (t166 * t180 + t170 * t174) * t167 + t166 * t184;
	t158 = t173 * t187 + (t173 * t182 + t174 * t175) * t168;
	t157 = -t175 * t187 + (t173 * t174 - t175 * t182) * t168;
	t156 = -t166 * t183 + t177 * t170;
	t155 = -t177 * t166 - t170 * t183;
	t154 = t156 * t175 - t163 * t173;
	t153 = t156 * t173 + t163 * t175;
	t152 = t155 * t175 + t173 * t162;
	t151 = t155 * t173 - t175 * t162;
	t1 = [t151 * t169 + t159 * t165, -t151 * t165 + t159 * t169, -t152, t151 * pkin(3) - t152 * qJ(4) + (t172 * t179 + t188) * t176 + (-t172 * t189 + t178) * t174 + t166 * t185 + t170 * pkin(1) + 0; t153 * t169 + t160 * t165, -t153 * t165 + t160 * t169, -t154, t153 * pkin(3) - t154 * qJ(4) + (-t172 * t178 + t189) * t176 + (t172 * t188 + t179) * t174 - t170 * t185 + t166 * pkin(1) + 0; t158 * t169 + t161 * t165, -t158 * t165 + t161 * t169, t157, t158 * pkin(3) + t157 * qJ(4) + t164 * t172 + qJ(1) + 0 + (pkin(2) * t174 - pkin(9) * t186) * t168; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:12:14
	% EndTime: 2020-11-04 21:12:14
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (117->55), mult. (262->89), div. (0->0), fcn. (342->14), ass. (0->45)
	t235 = sin(pkin(13)) * pkin(4);
	t211 = sin(pkin(7));
	t234 = pkin(9) * t211;
	t210 = sin(pkin(12));
	t233 = t210 * pkin(2);
	t213 = cos(pkin(12));
	t232 = t213 * pkin(2);
	t215 = cos(pkin(6));
	t231 = t211 * t215;
	t220 = cos(qJ(2));
	t230 = t211 * t220;
	t214 = cos(pkin(7));
	t205 = t214 * pkin(9) + pkin(8);
	t212 = sin(pkin(6));
	t229 = t212 * t205;
	t228 = t212 * t214;
	t218 = sin(qJ(2));
	t227 = t214 * t218;
	t226 = t214 * t220;
	t225 = t215 * t218;
	t224 = t215 * t220;
	t223 = t210 * t234;
	t222 = t213 * t234;
	t221 = -t211 * t212 + t214 * t224;
	t219 = cos(qJ(3));
	t217 = sin(qJ(3));
	t216 = -pkin(10) - qJ(4);
	t208 = pkin(13) + qJ(5);
	t207 = cos(t208);
	t206 = sin(t208);
	t204 = cos(pkin(13)) * pkin(4) + pkin(3);
	t203 = t210 * t220 + t213 * t225;
	t202 = t210 * t225 - t213 * t220;
	t201 = -t212 * t230 + t215 * t214;
	t200 = -(-t210 * t218 + t213 * t224) * t211 - t213 * t228;
	t199 = (t210 * t224 + t213 * t218) * t211 + t210 * t228;
	t198 = t217 * t231 + (t217 * t226 + t218 * t219) * t212;
	t197 = -t219 * t231 + (t217 * t218 - t219 * t226) * t212;
	t196 = -t210 * t227 + t221 * t213;
	t195 = -t221 * t210 - t213 * t227;
	t194 = t196 * t219 - t203 * t217;
	t193 = t196 * t217 + t203 * t219;
	t192 = t195 * t219 + t217 * t202;
	t191 = t195 * t217 - t219 * t202;
	t1 = [t191 * t207 + t199 * t206, -t191 * t206 + t199 * t207, -t192, t191 * t204 + t192 * t216 + t199 * t235 + (t215 * t223 + t232) * t220 + (-t215 * t233 + t222) * t218 + t210 * t229 + t213 * pkin(1) + 0; t193 * t207 + t200 * t206, -t193 * t206 + t200 * t207, -t194, t193 * t204 + t194 * t216 + t200 * t235 + (-t215 * t222 + t233) * t220 + (t215 * t232 + t223) * t218 - t213 * t229 + t210 * pkin(1) + 0; t198 * t207 + t201 * t206, -t198 * t206 + t201 * t207, t197, t201 * t235 - t197 * t216 + t198 * t204 + t205 * t215 + qJ(1) + 0 + (pkin(2) * t218 - pkin(9) * t230) * t212; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:12:14
	% EndTime: 2020-11-04 21:12:14
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (211->65), mult. (467->107), div. (0->0), fcn. (617->16), ass. (0->53)
	t289 = sin(pkin(13)) * pkin(4);
	t262 = sin(pkin(7));
	t288 = pkin(9) * t262;
	t261 = sin(pkin(12));
	t287 = t261 * pkin(2);
	t264 = cos(pkin(12));
	t286 = t264 * pkin(2);
	t266 = cos(pkin(6));
	t285 = t262 * t266;
	t273 = cos(qJ(2));
	t284 = t262 * t273;
	t265 = cos(pkin(7));
	t256 = t265 * pkin(9) + pkin(8);
	t263 = sin(pkin(6));
	t283 = t263 * t256;
	t282 = t263 * t265;
	t270 = sin(qJ(2));
	t281 = t265 * t270;
	t280 = t265 * t273;
	t279 = t266 * t270;
	t278 = t266 * t273;
	t277 = t261 * t288;
	t276 = t264 * t288;
	t275 = t261 * t273 + t264 * t279;
	t274 = -t262 * t263 + t265 * t278;
	t272 = cos(qJ(3));
	t271 = cos(qJ(6));
	t269 = sin(qJ(3));
	t268 = sin(qJ(6));
	t267 = -pkin(10) - qJ(4);
	t259 = pkin(13) + qJ(5);
	t258 = cos(t259);
	t257 = sin(t259);
	t255 = cos(pkin(13)) * pkin(4) + pkin(3);
	t254 = t261 * t279 - t264 * t273;
	t253 = -t263 * t284 + t266 * t265;
	t251 = -(-t261 * t270 + t264 * t278) * t262 - t264 * t282;
	t250 = (t261 * t278 + t264 * t270) * t262 + t261 * t282;
	t249 = t269 * t285 + (t269 * t280 + t270 * t272) * t263;
	t248 = -t272 * t285 + (t269 * t270 - t272 * t280) * t263;
	t247 = -t261 * t281 + t274 * t264;
	t246 = -t274 * t261 - t264 * t281;
	t245 = t247 * t272 - t275 * t269;
	t244 = t247 * t269 + t275 * t272;
	t243 = t246 * t272 + t269 * t254;
	t242 = t246 * t269 - t272 * t254;
	t241 = t249 * t258 + t253 * t257;
	t240 = t249 * t257 - t253 * t258;
	t239 = t244 * t258 + t251 * t257;
	t238 = t244 * t257 - t251 * t258;
	t237 = t242 * t258 + t250 * t257;
	t236 = t242 * t257 - t250 * t258;
	t1 = [t237 * t271 - t243 * t268, -t237 * t268 - t243 * t271, t236, t237 * pkin(5) + t236 * pkin(11) + t242 * t255 + t243 * t267 + t250 * t289 + (t266 * t277 + t286) * t273 + (-t266 * t287 + t276) * t270 + t261 * t283 + t264 * pkin(1) + 0; t239 * t271 - t245 * t268, -t239 * t268 - t245 * t271, t238, t239 * pkin(5) + t238 * pkin(11) + t244 * t255 + t245 * t267 + t251 * t289 + (-t266 * t276 + t287) * t273 + (t266 * t286 + t277) * t270 - t264 * t283 + t261 * pkin(1) + 0; t241 * t271 + t248 * t268, -t241 * t268 + t248 * t271, t240, t253 * t289 + t241 * pkin(5) + t240 * pkin(11) - t248 * t267 + t249 * t255 + t256 * t266 + qJ(1) + 0 + (pkin(2) * t270 - pkin(9) * t284) * t263; 0, 0, 0, 1;];
	Tc_mdh = t1;
end