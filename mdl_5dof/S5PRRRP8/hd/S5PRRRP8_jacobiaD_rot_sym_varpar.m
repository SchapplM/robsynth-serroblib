% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5PRRRP8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRRP8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:01:48
	% EndTime: 2019-12-05 17:01:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:01:48
	% EndTime: 2019-12-05 17:01:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:01:48
	% EndTime: 2019-12-05 17:01:48
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(9));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(9)) * cos(pkin(5));
	t37 = t39 * t42 - t41 * t45;
	t34 = 0.1e1 / t37 ^ 2;
	t49 = qJD(2) * t34;
	t36 = t39 * t41 + t42 * t45;
	t33 = t36 ^ 2;
	t30 = t33 * t34 + 0.1e1;
	t46 = t37 * t49;
	t47 = t36 / t37 * t49;
	t48 = (t33 * t47 + t36 * t46) / t30 ^ 2;
	t28 = 0.1e1 / t30;
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, -0.2e1 * t48 + 0.2e1 * (t28 * t46 + (t28 * t47 - t34 * t48) * t36) * t36, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:01:49
	% EndTime: 2019-12-05 17:01:49
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (756->55), mult. (2271->133), div. (423->14), fcn. (2956->11), ass. (0->65)
	t139 = sin(qJ(2));
	t141 = cos(qJ(2));
	t136 = sin(pkin(9));
	t166 = cos(pkin(5));
	t155 = t136 * t166;
	t165 = cos(pkin(9));
	t127 = -t139 * t155 + t165 * t141;
	t138 = sin(qJ(3));
	t140 = cos(qJ(3));
	t137 = sin(pkin(5));
	t160 = t136 * t137;
	t149 = -t127 * t138 + t140 * t160;
	t170 = t149 * qJD(3);
	t152 = t166 * t165;
	t123 = t136 * t139 - t141 * t152;
	t159 = t137 * t141;
	t113 = atan2(-t123, -t159);
	t111 = sin(t113);
	t112 = cos(t113);
	t98 = -t111 * t123 - t112 * t159;
	t95 = 0.1e1 / t98;
	t110 = t127 * t140 + t138 * t160;
	t106 = 0.1e1 / t110;
	t133 = 0.1e1 / t141;
	t107 = 0.1e1 / t110 ^ 2;
	t134 = 0.1e1 / t141 ^ 2;
	t96 = 0.1e1 / t98 ^ 2;
	t105 = t149 ^ 2;
	t102 = t105 * t107 + 0.1e1;
	t148 = -t165 * t139 - t141 * t155;
	t119 = t148 * qJD(2);
	t103 = t110 * qJD(3) + t119 * t138;
	t163 = t107 * t149;
	t104 = t119 * t140 + t170;
	t164 = t104 * t106 * t107;
	t169 = 0.1e1 / t102 ^ 2 * (-t103 * t163 - t105 * t164);
	t125 = t136 * t141 + t139 * t152;
	t161 = t134 * t139;
	t156 = t123 * t161;
	t150 = t125 * t133 + t156;
	t121 = t123 ^ 2;
	t132 = 0.1e1 / t137 ^ 2;
	t116 = t121 * t132 * t134 + 0.1e1;
	t114 = 0.1e1 / t116;
	t131 = 0.1e1 / t137;
	t162 = t114 * t131;
	t91 = t150 * t162;
	t168 = t123 * t91;
	t167 = t148 * t96;
	t158 = qJD(2) * t139;
	t157 = -0.2e1 * t169;
	t151 = -t106 * t138 - t140 * t163;
	t135 = t133 * t134;
	t122 = t148 ^ 2;
	t120 = t127 * qJD(2);
	t118 = t125 * qJD(2);
	t117 = t123 * qJD(2);
	t100 = 0.1e1 / t102;
	t97 = t95 * t96;
	t94 = t122 * t96 + 0.1e1;
	t90 = (qJD(2) * t156 + t118 * t133) * t162;
	t88 = (t137 * t139 - t168) * t112 + (t91 * t159 - t125) * t111;
	t87 = (-t123 * t90 + t137 * t158) * t112 + (t90 * t159 - t118) * t111;
	t86 = (-0.2e1 * t150 * (t118 * t123 * t134 + t121 * t135 * t158) * t132 / t116 ^ 2 + (t118 * t161 - t117 * t133 + (t125 * t161 + (0.2e1 * t135 * t139 ^ 2 + t133) * t123) * qJD(2)) * t114) * t131;
	t1 = [0, t86, 0, 0, 0; 0, 0.2e1 * (-t127 * t95 - t88 * t167) / t94 ^ 2 * (-t122 * t87 * t97 - t120 * t167) + (-t88 * t120 * t96 + t119 * t95 + (-0.2e1 * t148 * t88 * t97 - t127 * t96) * t87 - (-(-t118 * t91 - t123 * t86 - t125 * t90 + (t90 * t91 + qJD(2)) * t159) * t112 - (t90 * t168 + t117 + (t141 * t86 + (-qJD(2) * t91 - t90) * t139) * t137) * t111) * t167) / t94, 0, 0, 0; 0, -t151 * t148 * t157 + (t151 * t120 - ((-qJD(3) * t106 + 0.2e1 * t149 * t164) * t140 + (t103 * t140 + (t104 + t170) * t138) * t107) * t148) * t100, t157 - 0.2e1 * (t100 * t103 * t107 - (-t100 * t164 - t107 * t169) * t149) * t149, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:01:49
	% EndTime: 2019-12-05 17:01:50
	% DurationCPUTime: 0.98s
	% Computational Cost: add. (3002->109), mult. (9085->226), div. (559->12), fcn. (11668->13), ass. (0->103)
	t207 = sin(pkin(9));
	t209 = cos(pkin(9));
	t213 = sin(qJ(2));
	t210 = cos(pkin(5));
	t216 = cos(qJ(2));
	t242 = t210 * t216;
	t197 = -t207 * t213 + t209 * t242;
	t190 = t197 * qJD(2);
	t243 = t210 * t213;
	t198 = t207 * t216 + t209 * t243;
	t212 = sin(qJ(3));
	t208 = sin(pkin(5));
	t246 = t208 * t212;
	t232 = t209 * t246;
	t215 = cos(qJ(3));
	t239 = qJD(3) * t215;
	t162 = -qJD(3) * t232 + t190 * t212 + t198 * t239;
	t245 = t208 * t215;
	t182 = t198 * t212 + t209 * t245;
	t180 = t182 ^ 2;
	t201 = -t210 * t215 + t213 * t246;
	t195 = 0.1e1 / t201 ^ 2;
	t176 = t180 * t195 + 0.1e1;
	t174 = 0.1e1 / t176;
	t202 = t210 * t212 + t213 * t245;
	t240 = qJD(2) * t216;
	t231 = t208 * t240;
	t187 = t202 * qJD(3) + t212 * t231;
	t194 = 0.1e1 / t201;
	t250 = t182 * t195;
	t146 = (-t162 * t194 + t187 * t250) * t174;
	t177 = atan2(-t182, t201);
	t172 = sin(t177);
	t173 = cos(t177);
	t228 = -t172 * t201 - t173 * t182;
	t142 = t228 * t146 - t172 * t162 + t173 * t187;
	t156 = -t172 * t182 + t173 * t201;
	t153 = 0.1e1 / t156;
	t154 = 0.1e1 / t156 ^ 2;
	t263 = t142 * t153 * t154;
	t233 = t207 * t243;
	t200 = t209 * t216 - t233;
	t225 = -t200 * t212 + t207 * t245;
	t262 = -0.2e1 * t225 * t263;
	t244 = t208 * t216;
	t224 = -t194 * t197 + t244 * t250;
	t261 = t212 * t224;
	t249 = t187 * t194 * t195;
	t260 = -0.2e1 * (t162 * t250 - t180 * t249) / t176 ^ 2;
	t186 = t200 * t215 + t207 * t246;
	t199 = t207 * t242 + t209 * t213;
	t211 = sin(qJ(4));
	t214 = cos(qJ(4));
	t171 = t186 * t214 + t199 * t211;
	t167 = 0.1e1 / t171;
	t168 = 0.1e1 / t171 ^ 2;
	t192 = t199 * qJD(2);
	t165 = t225 * qJD(3) - t192 * t215;
	t193 = -qJD(2) * t233 + t209 * t240;
	t157 = t171 * qJD(4) + t165 * t211 - t193 * t214;
	t170 = t186 * t211 - t199 * t214;
	t166 = t170 ^ 2;
	t161 = t166 * t168 + 0.1e1;
	t254 = t168 * t170;
	t238 = qJD(4) * t170;
	t158 = t165 * t214 + t193 * t211 - t238;
	t257 = t158 * t167 * t168;
	t259 = (t157 * t254 - t166 * t257) / t161 ^ 2;
	t258 = t154 * t225;
	t164 = t186 * qJD(3) - t192 * t212;
	t256 = t164 * t154;
	t255 = t167 * t211;
	t253 = t170 * t214;
	t252 = t172 * t225;
	t251 = t173 * t225;
	t248 = t199 * t212;
	t247 = t199 * t215;
	t241 = qJD(2) * t213;
	t181 = t225 ^ 2;
	t152 = t181 * t154 + 0.1e1;
	t237 = 0.2e1 * (-t181 * t263 - t225 * t256) / t152 ^ 2;
	t236 = -0.2e1 * t259;
	t234 = t170 * t257;
	t230 = -0.2e1 * t182 * t249;
	t229 = qJD(4) * t247 - t192;
	t227 = t168 * t253 - t255;
	t184 = t198 * t215 - t232;
	t226 = -t184 * t194 + t202 * t250;
	t223 = qJD(3) * t248 + qJD(4) * t200 - t193 * t215;
	t191 = t198 * qJD(2);
	t188 = -t201 * qJD(3) + t215 * t231;
	t179 = t200 * t211 - t214 * t247;
	t178 = -t200 * t214 - t211 * t247;
	t163 = -t182 * qJD(3) + t190 * t215;
	t159 = 0.1e1 / t161;
	t149 = 0.1e1 / t152;
	t148 = t174 * t261;
	t147 = t226 * t174;
	t144 = (-t172 * t197 + t173 * t244) * t212 + t228 * t148;
	t143 = t228 * t147 - t172 * t184 + t173 * t202;
	t141 = t226 * t260 + (t202 * t230 - t163 * t194 + (t162 * t202 + t182 * t188 + t184 * t187) * t195) * t174;
	t139 = t260 * t261 + (t224 * t239 + (t230 * t244 + t191 * t194 + (t187 * t197 + (t162 * t216 - t182 * t241) * t208) * t195) * t212) * t174;
	t1 = [0, t139, t141, 0, 0; 0, (-t144 * t258 + t153 * t248) * t237 + ((-t193 * t212 - t199 * t239) * t153 + (-t256 + t262) * t144 + (t248 * t142 + (-t139 * t182 - t148 * t162 + (-t212 * t241 + t216 * t239) * t208 + (-t148 * t201 - t197 * t212) * t146) * t251 + (-t197 * t239 - t139 * t201 - t148 * t187 + t191 * t212 + (t148 * t182 - t212 * t244) * t146) * t252) * t154) * t149, (-t143 * t258 - t153 * t186) * t237 + (t143 * t262 + t165 * t153 + (-t186 * t142 - t143 * t164 + (-t141 * t182 - t147 * t162 + t188 + (-t147 * t201 - t184) * t146) * t251 + (-t141 * t201 - t147 * t187 - t163 + (t147 * t182 - t202) * t146) * t252) * t154) * t149, 0, 0; 0, 0.2e1 * (-t167 * t178 + t179 * t254) * t259 + (0.2e1 * t179 * t234 - t229 * t167 * t214 + t223 * t255 + (-t229 * t170 * t211 - t179 * t157 - t178 * t158 - t223 * t253) * t168) * t159, -t227 * t225 * t236 + (t227 * t164 - ((-qJD(4) * t167 - 0.2e1 * t234) * t214 + (t157 * t214 + (t158 - t238) * t211) * t168) * t225) * t159, t236 + 0.2e1 * (t157 * t168 * t159 + (-t159 * t257 - t168 * t259) * t170) * t170, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:01:49
	% EndTime: 2019-12-05 17:01:50
	% DurationCPUTime: 1.74s
	% Computational Cost: add. (7242->156), mult. (21125->307), div. (784->12), fcn. (27225->13), ass. (0->123)
	t251 = sin(qJ(3));
	t254 = cos(qJ(3));
	t252 = sin(qJ(2));
	t255 = cos(qJ(2));
	t312 = cos(pkin(9));
	t313 = cos(pkin(5));
	t276 = t313 * t312;
	t311 = sin(pkin(9));
	t264 = -t252 * t276 - t311 * t255;
	t249 = sin(pkin(5));
	t281 = t249 * t312;
	t223 = t251 * t264 - t254 * t281;
	t274 = t311 * t252 - t255 * t276;
	t239 = t274 * qJD(2);
	t206 = t223 * qJD(3) - t239 * t254;
	t250 = sin(qJ(4));
	t253 = cos(qJ(4));
	t266 = t251 * t281 + t254 * t264;
	t216 = t274 * t250 - t253 * t266;
	t262 = t264 * qJD(2);
	t188 = t216 * qJD(4) + t206 * t250 + t253 * t262;
	t267 = t274 * t253;
	t214 = -t250 * t266 - t267;
	t209 = t214 ^ 2;
	t298 = t249 * t252;
	t246 = t313 * t251 + t254 * t298;
	t295 = t253 * t255;
	t232 = t246 * t250 + t249 * t295;
	t228 = 0.1e1 / t232 ^ 2;
	t200 = t209 * t228 + 0.1e1;
	t198 = 0.1e1 / t200;
	t245 = -t251 * t298 + t313 * t254;
	t294 = qJD(2) * t249;
	t282 = t255 * t294;
	t231 = t245 * qJD(3) + t254 * t282;
	t297 = t250 * t255;
	t233 = t246 * t253 - t249 * t297;
	t283 = t252 * t294;
	t202 = t233 * qJD(4) + t231 * t250 - t253 * t283;
	t227 = 0.1e1 / t232;
	t302 = t214 * t228;
	t175 = (-t188 * t227 + t202 * t302) * t198;
	t201 = atan2(-t214, t232);
	t195 = sin(t201);
	t196 = cos(t201);
	t273 = -t195 * t232 - t196 * t214;
	t171 = t273 * t175 - t195 * t188 + t196 * t202;
	t187 = -t195 * t214 + t196 * t232;
	t184 = 0.1e1 / t187;
	t185 = 0.1e1 / t187 ^ 2;
	t317 = t171 * t184 * t185;
	t275 = t313 * t311;
	t263 = -t312 * t252 - t255 * t275;
	t244 = -t252 * t275 + t312 * t255;
	t280 = t249 * t311;
	t265 = -t244 * t254 - t251 * t280;
	t217 = -t250 * t265 + t253 * t263;
	t316 = -0.2e1 * t217;
	t278 = 0.2e1 * t217 * t317;
	t269 = -t223 * t227 + t245 * t302;
	t315 = t250 * t269;
	t304 = t202 * t227 * t228;
	t314 = -0.2e1 * (t188 * t302 - t209 * t304) / t200 ^ 2;
	t218 = -t250 * t263 - t253 * t265;
	t211 = 0.1e1 / t218;
	t212 = 0.1e1 / t218 ^ 2;
	t225 = -t244 * t251 + t254 * t280;
	t222 = t225 ^ 2;
	t301 = t222 * t212;
	t197 = 0.1e1 + t301;
	t240 = t263 * qJD(2);
	t207 = t265 * qJD(3) - t240 * t251;
	t208 = t225 * qJD(3) + t240 * t254;
	t241 = t244 * qJD(2);
	t191 = -t217 * qJD(4) + t208 * t253 + t241 * t250;
	t307 = t191 * t211 * t212;
	t284 = t222 * t307;
	t303 = t212 * t225;
	t310 = (t207 * t303 - t284) / t197 ^ 2;
	t309 = t185 * t217;
	t190 = t218 * qJD(4) + t208 * t250 - t241 * t253;
	t308 = t190 * t185;
	t306 = t195 * t217;
	t305 = t196 * t217;
	t300 = t225 * t250;
	t299 = t263 * t254;
	t296 = t253 * t254;
	t293 = qJD(2) * t254;
	t292 = qJD(3) * t251;
	t291 = qJD(4) * t250;
	t290 = qJD(4) * t253;
	t289 = t254 * qJD(4);
	t210 = t217 ^ 2;
	t183 = t210 * t185 + 0.1e1;
	t288 = 0.2e1 * (-t210 * t317 + t217 * t308) / t183 ^ 2;
	t287 = 0.2e1 * t310;
	t285 = t225 * t307;
	t277 = -0.2e1 * t214 * t304;
	t271 = -t216 * t227 + t233 * t302;
	t268 = t254 * t274;
	t219 = -t250 * t268 + t253 * t264;
	t235 = (-t252 * t253 + t254 * t297) * t249;
	t270 = -t219 * t227 + t235 * t302;
	t230 = -t246 * qJD(3) - t251 * t282;
	t221 = t244 * t250 + t263 * t296;
	t220 = -t244 * t253 + t250 * t299;
	t205 = t266 * qJD(3) + t239 * t251;
	t204 = ((-qJD(2) + t289) * t295 + (-t255 * t292 + (qJD(4) - t293) * t252) * t250) * t249;
	t203 = -t232 * qJD(4) + t231 * t253 + t250 * t283;
	t193 = 0.1e1 / t197;
	t192 = t239 * t253 - t264 * t291 - t268 * t290 + (t264 * t293 + t274 * t292) * t250;
	t189 = qJD(4) * t267 + t206 * t253 - t250 * t262 + t266 * t291;
	t181 = 0.1e1 / t183;
	t180 = t198 * t315;
	t179 = t270 * t198;
	t177 = t271 * t198;
	t174 = (-t195 * t223 + t196 * t245) * t250 + t273 * t180;
	t173 = t273 * t179 - t195 * t219 + t196 * t235;
	t172 = t273 * t177 - t195 * t216 + t196 * t233;
	t170 = t270 * t314 + (t235 * t277 - t192 * t227 + (t188 * t235 + t202 * t219 + t204 * t214) * t228) * t198;
	t168 = t271 * t314 + (t233 * t277 - t189 * t227 + (t188 * t233 + t202 * t216 + t203 * t214) * t228) * t198;
	t167 = t314 * t315 + (t269 * t290 + (t245 * t277 - t205 * t227 + (t188 * t245 + t202 * t223 + t214 * t230) * t228) * t250) * t198;
	t1 = [0, t170, t167, t168, 0; 0, (t173 * t309 - t184 * t220) * t288 + (t173 * t278 + (-t220 * t171 - t173 * t190 - (-t170 * t214 - t179 * t188 + t204 + (-t179 * t232 - t219) * t175) * t305 - (-t170 * t232 - t179 * t202 - t192 + (t179 * t214 - t235) * t175) * t306) * t185 + ((t263 * t289 - t240) * t253 + (qJD(4) * t244 - t241 * t254 - t263 * t292) * t250) * t184) * t181, (t174 * t309 - t184 * t300) * t288 + ((t207 * t250 + t225 * t290) * t184 + (-t308 + t278) * t174 + (-t300 * t171 - (t245 * t290 - t167 * t214 - t180 * t188 + t230 * t250 + (-t180 * t232 - t223 * t250) * t175) * t305 - (-t223 * t290 - t167 * t232 - t180 * t202 - t205 * t250 + (t180 * t214 - t245 * t250) * t175) * t306) * t185) * t181, (t172 * t309 - t184 * t218) * t288 + (t172 * t278 + t191 * t184 + (-t218 * t171 - t172 * t190 - (-t168 * t214 - t177 * t188 + t203 + (-t177 * t232 - t216) * t175) * t305 - (-t168 * t232 - t177 * t202 - t189 + (t177 * t214 - t233) * t175) * t306) * t185) * t181, 0; 0, (t211 * t251 * t263 + t221 * t303) * t287 + (0.2e1 * t221 * t285 + (-qJD(3) * t299 + t241 * t251) * t211 + (-(t240 * t250 - t241 * t296 + t244 * t290) * t225 - t221 * t207 - (-t251 * t191 - (t250 * t289 + t253 * t292) * t225) * t263) * t212) * t193, (-t211 * t265 + t253 * t301) * t287 + (0.2e1 * t253 * t284 - t208 * t211 + (-0.2e1 * t207 * t225 * t253 - t191 * t265 + t222 * t291) * t212) * t193, t303 * t310 * t316 + (t285 * t316 + (t190 * t225 + t207 * t217) * t212) * t193, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end