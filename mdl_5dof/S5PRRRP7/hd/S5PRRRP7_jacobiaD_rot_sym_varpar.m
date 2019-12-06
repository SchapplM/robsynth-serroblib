% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRP7
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
%   Wie in S5PRRRP7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRRP7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:57:06
	% EndTime: 2019-12-05 16:57:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:57:06
	% EndTime: 2019-12-05 16:57:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:57:06
	% EndTime: 2019-12-05 16:57:07
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
	% StartTime: 2019-12-05 16:57:06
	% EndTime: 2019-12-05 16:57:07
	% DurationCPUTime: 0.41s
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
	% StartTime: 2019-12-05 16:57:07
	% EndTime: 2019-12-05 16:57:08
	% DurationCPUTime: 0.97s
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
	% StartTime: 2019-12-05 16:57:07
	% EndTime: 2019-12-05 16:57:08
	% DurationCPUTime: 0.95s
	% Computational Cost: add. (3002->109), mult. (9085->226), div. (559->12), fcn. (11668->13), ass. (0->102)
	t215 = sin(pkin(9));
	t217 = cos(pkin(9));
	t221 = sin(qJ(2));
	t218 = cos(pkin(5));
	t224 = cos(qJ(2));
	t250 = t218 * t224;
	t205 = -t215 * t221 + t217 * t250;
	t198 = t205 * qJD(2);
	t251 = t218 * t221;
	t206 = t215 * t224 + t217 * t251;
	t220 = sin(qJ(3));
	t216 = sin(pkin(5));
	t254 = t216 * t220;
	t240 = t217 * t254;
	t223 = cos(qJ(3));
	t247 = qJD(3) * t223;
	t170 = -qJD(3) * t240 + t198 * t220 + t206 * t247;
	t253 = t216 * t223;
	t190 = t206 * t220 + t217 * t253;
	t188 = t190 ^ 2;
	t209 = -t218 * t223 + t221 * t254;
	t203 = 0.1e1 / t209 ^ 2;
	t184 = t188 * t203 + 0.1e1;
	t182 = 0.1e1 / t184;
	t210 = t218 * t220 + t221 * t253;
	t248 = qJD(2) * t224;
	t239 = t216 * t248;
	t195 = t210 * qJD(3) + t220 * t239;
	t202 = 0.1e1 / t209;
	t258 = t190 * t203;
	t154 = (-t170 * t202 + t195 * t258) * t182;
	t185 = atan2(-t190, t209);
	t180 = sin(t185);
	t181 = cos(t185);
	t236 = -t180 * t209 - t181 * t190;
	t150 = t236 * t154 - t180 * t170 + t181 * t195;
	t164 = -t180 * t190 + t181 * t209;
	t161 = 0.1e1 / t164;
	t162 = 0.1e1 / t164 ^ 2;
	t270 = t150 * t161 * t162;
	t241 = t215 * t251;
	t208 = t217 * t224 - t241;
	t233 = -t208 * t220 + t215 * t253;
	t269 = -0.2e1 * t233 * t270;
	t252 = t216 * t224;
	t232 = -t202 * t205 + t252 * t258;
	t268 = t220 * t232;
	t257 = t195 * t202 * t203;
	t267 = -0.2e1 * (t170 * t258 - t188 * t257) / t184 ^ 2;
	t194 = t208 * t223 + t215 * t254;
	t207 = t215 * t250 + t217 * t221;
	t219 = sin(qJ(4));
	t222 = cos(qJ(4));
	t179 = t194 * t222 + t207 * t219;
	t175 = 0.1e1 / t179;
	t176 = 0.1e1 / t179 ^ 2;
	t200 = t207 * qJD(2);
	t173 = t233 * qJD(3) - t200 * t223;
	t201 = -qJD(2) * t241 + t217 * t248;
	t165 = t179 * qJD(4) + t173 * t219 - t201 * t222;
	t178 = t194 * t219 - t207 * t222;
	t174 = t178 ^ 2;
	t169 = t174 * t176 + 0.1e1;
	t262 = t176 * t178;
	t246 = qJD(4) * t178;
	t166 = t173 * t222 + t201 * t219 - t246;
	t264 = t166 * t175 * t176;
	t266 = (t165 * t262 - t174 * t264) / t169 ^ 2;
	t265 = t162 * t233;
	t263 = t175 * t219;
	t261 = t178 * t222;
	t260 = t180 * t233;
	t259 = t181 * t233;
	t256 = t207 * t220;
	t255 = t207 * t223;
	t249 = qJD(2) * t221;
	t189 = t233 ^ 2;
	t160 = t189 * t162 + 0.1e1;
	t172 = t194 * qJD(3) - t200 * t220;
	t245 = 0.2e1 * (-t172 * t265 - t189 * t270) / t160 ^ 2;
	t244 = -0.2e1 * t266;
	t242 = t178 * t264;
	t238 = -0.2e1 * t190 * t257;
	t237 = qJD(4) * t255 - t200;
	t235 = t176 * t261 - t263;
	t192 = t206 * t223 - t240;
	t234 = -t192 * t202 + t210 * t258;
	t231 = qJD(3) * t256 + qJD(4) * t208 - t201 * t223;
	t199 = t206 * qJD(2);
	t196 = -t209 * qJD(3) + t223 * t239;
	t187 = t208 * t219 - t222 * t255;
	t186 = -t208 * t222 - t219 * t255;
	t171 = -t190 * qJD(3) + t198 * t223;
	t167 = 0.1e1 / t169;
	t157 = 0.1e1 / t160;
	t156 = t182 * t268;
	t155 = t234 * t182;
	t152 = (-t180 * t205 + t181 * t252) * t220 + t236 * t156;
	t151 = t236 * t155 - t180 * t192 + t181 * t210;
	t149 = t234 * t267 + (t210 * t238 - t171 * t202 + (t170 * t210 + t190 * t196 + t192 * t195) * t203) * t182;
	t147 = t267 * t268 + (t232 * t247 + (t238 * t252 + t199 * t202 + (t195 * t205 + (t170 * t224 - t190 * t249) * t216) * t203) * t220) * t182;
	t1 = [0, t147, t149, 0, 0; 0, (-t152 * t265 + t161 * t256) * t245 + ((-t201 * t220 - t207 * t247) * t161 + t152 * t269 + (-t152 * t172 + t256 * t150 + (-t147 * t190 - t156 * t170 + (-t220 * t249 + t224 * t247) * t216 + (-t156 * t209 - t205 * t220) * t154) * t259 + (-t205 * t247 - t147 * t209 - t156 * t195 + t199 * t220 + (t156 * t190 - t220 * t252) * t154) * t260) * t162) * t157, (-t151 * t265 - t161 * t194) * t245 + (t151 * t269 + t173 * t161 + (-t194 * t150 - t151 * t172 + (-t149 * t190 - t155 * t170 + t196 + (-t155 * t209 - t192) * t154) * t259 + (-t149 * t209 - t155 * t195 - t171 + (t155 * t190 - t210) * t154) * t260) * t162) * t157, 0, 0; 0, 0.2e1 * (-t175 * t186 + t187 * t262) * t266 + (0.2e1 * t187 * t242 - t237 * t175 * t222 + t231 * t263 + (-t237 * t178 * t219 - t187 * t165 - t186 * t166 - t231 * t261) * t176) * t167, -t235 * t233 * t244 + (t235 * t172 - ((-qJD(4) * t175 - 0.2e1 * t242) * t222 + (t165 * t222 + (t166 - t246) * t219) * t176) * t233) * t167, t244 + 0.2e1 * (t165 * t176 * t167 + (-t167 * t264 - t176 * t266) * t178) * t178, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end