% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6PRRPPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPPR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(10));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(10)) * cos(pkin(6));
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
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, -0.2e1 * t48 + 0.2e1 * (t28 * t46 + (t28 * t47 - t34 * t48) * t36) * t36, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.63s
	% Computational Cost: add. (756->55), mult. (2271->133), div. (423->14), fcn. (2956->11), ass. (0->65)
	t139 = sin(qJ(2));
	t141 = cos(qJ(2));
	t136 = sin(pkin(10));
	t166 = cos(pkin(6));
	t155 = t136 * t166;
	t165 = cos(pkin(10));
	t127 = -t139 * t155 + t165 * t141;
	t138 = sin(qJ(3));
	t140 = cos(qJ(3));
	t137 = sin(pkin(6));
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
	t117 = qJD(2) * t123;
	t100 = 0.1e1 / t102;
	t97 = t95 * t96;
	t94 = t122 * t96 + 0.1e1;
	t90 = (qJD(2) * t156 + t118 * t133) * t162;
	t88 = (t137 * t139 - t168) * t112 + (t91 * t159 - t125) * t111;
	t87 = (-t123 * t90 + t137 * t158) * t112 + (t90 * t159 - t118) * t111;
	t86 = (-0.2e1 * t150 * (t118 * t123 * t134 + t121 * t135 * t158) * t132 / t116 ^ 2 + (t118 * t161 - t117 * t133 + (t125 * t161 + (0.2e1 * t135 * t139 ^ 2 + t133) * t123) * qJD(2)) * t114) * t131;
	t1 = [0, t86, 0, 0, 0, 0; 0, 0.2e1 * (-t127 * t95 - t88 * t167) / t94 ^ 2 * (-t122 * t97 * t87 - t120 * t167) + (-t88 * t120 * t96 + t119 * t95 + (-0.2e1 * t148 * t88 * t97 - t127 * t96) * t87 - (-(-t118 * t91 - t123 * t86 - t125 * t90 + (t90 * t91 + qJD(2)) * t159) * t112 - (t90 * t168 + t117 + (t141 * t86 + (-qJD(2) * t91 - t90) * t139) * t137) * t111) * t167) / t94, 0, 0, 0, 0; 0, -t151 * t148 * t157 + (t151 * t120 - ((-qJD(3) * t106 + 0.2e1 * t149 * t164) * t140 + (t103 * t140 + (t104 + t170) * t138) * t107) * t148) * t100, t157 - 0.2e1 * (t100 * t103 * t107 - (-t100 * t164 - t107 * t169) * t149) * t149, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:45
	% DurationCPUTime: 1.37s
	% Computational Cost: add. (2688->101), mult. (8196->213), div. (535->12), fcn. (10572->13), ass. (0->98)
	t187 = sin(pkin(10));
	t190 = cos(pkin(10));
	t193 = sin(qJ(2));
	t191 = cos(pkin(6));
	t195 = cos(qJ(2));
	t218 = t191 * t195;
	t176 = -t187 * t193 + t190 * t218;
	t169 = t176 * qJD(2);
	t219 = t191 * t193;
	t177 = t187 * t195 + t190 * t219;
	t192 = sin(qJ(3));
	t188 = sin(pkin(6));
	t222 = t188 * t192;
	t209 = t190 * t222;
	t194 = cos(qJ(3));
	t215 = qJD(3) * t194;
	t142 = -qJD(3) * t209 + t169 * t192 + t177 * t215;
	t221 = t188 * t194;
	t162 = t177 * t192 + t190 * t221;
	t160 = t162 ^ 2;
	t180 = -t191 * t194 + t193 * t222;
	t174 = 0.1e1 / t180 ^ 2;
	t156 = t160 * t174 + 0.1e1;
	t154 = 0.1e1 / t156;
	t181 = t191 * t192 + t193 * t221;
	t216 = qJD(2) * t195;
	t208 = t188 * t216;
	t167 = t181 * qJD(3) + t192 * t208;
	t173 = 0.1e1 / t180;
	t226 = t162 * t174;
	t126 = (-t142 * t173 + t167 * t226) * t154;
	t157 = atan2(-t162, t180);
	t152 = sin(t157);
	t153 = cos(t157);
	t206 = -t152 * t180 - t153 * t162;
	t122 = t206 * t126 - t142 * t152 + t153 * t167;
	t136 = -t152 * t162 + t153 * t180;
	t133 = 0.1e1 / t136;
	t134 = 0.1e1 / t136 ^ 2;
	t238 = t122 * t133 * t134;
	t210 = t187 * t219;
	t179 = t190 * t195 - t210;
	t166 = t179 * t194 + t187 * t222;
	t178 = t187 * t218 + t190 * t193;
	t186 = sin(pkin(11));
	t189 = cos(pkin(11));
	t150 = t166 * t186 - t178 * t189;
	t171 = t178 * qJD(2);
	t204 = -t179 * t192 + t187 * t221;
	t145 = t204 * qJD(3) - t171 * t194;
	t172 = -qJD(2) * t210 + t190 * t216;
	t141 = t145 * t189 + t172 * t186;
	t151 = t166 * t189 + t178 * t186;
	t147 = 0.1e1 / t151;
	t148 = 0.1e1 / t151 ^ 2;
	t232 = t141 * t147 * t148;
	t237 = 0.2e1 * t150 * t232;
	t236 = -0.2e1 * t204 * t238;
	t220 = t188 * t195;
	t202 = -t173 * t176 + t220 * t226;
	t235 = t192 * t202;
	t225 = t167 * t173 * t174;
	t234 = -0.2e1 * (t142 * t226 - t160 * t225) / t156 ^ 2;
	t233 = t134 * t204;
	t231 = t147 * t186;
	t230 = t148 * t150;
	t229 = t150 * t189;
	t228 = t152 * t204;
	t227 = t153 * t204;
	t224 = t178 * t192;
	t223 = t178 * t194;
	t217 = qJD(2) * t193;
	t161 = t204 ^ 2;
	t132 = t134 * t161 + 0.1e1;
	t144 = t166 * qJD(3) - t171 * t192;
	t214 = 0.2e1 * (-t144 * t233 - t161 * t238) / t132 ^ 2;
	t146 = t150 ^ 2;
	t139 = t146 * t148 + 0.1e1;
	t140 = t145 * t186 - t172 * t189;
	t213 = 0.2e1 * (t140 * t230 - t146 * t232) / t139 ^ 2;
	t207 = -0.2e1 * t162 * t225;
	t164 = t177 * t194 - t209;
	t205 = -t164 * t173 + t181 * t226;
	t203 = qJD(3) * t224 - t172 * t194;
	t170 = t177 * qJD(2);
	t168 = -t180 * qJD(3) + t194 * t208;
	t159 = t179 * t186 - t189 * t223;
	t158 = -t179 * t189 - t186 * t223;
	t143 = -t162 * qJD(3) + t169 * t194;
	t137 = 0.1e1 / t139;
	t129 = 0.1e1 / t132;
	t128 = t154 * t235;
	t127 = t205 * t154;
	t124 = (-t152 * t176 + t153 * t220) * t192 + t206 * t128;
	t123 = t206 * t127 - t152 * t164 + t153 * t181;
	t121 = t205 * t234 + (t181 * t207 - t143 * t173 + (t142 * t181 + t162 * t168 + t164 * t167) * t174) * t154;
	t119 = t234 * t235 + (t202 * t215 + (t207 * t220 + t170 * t173 + (t167 * t176 + (t142 * t195 - t162 * t217) * t188) * t174) * t192) * t154;
	t1 = [0, t119, t121, 0, 0, 0; 0, (-t124 * t233 + t133 * t224) * t214 + ((-t172 * t192 - t178 * t215) * t133 + t124 * t236 + (-t124 * t144 + t224 * t122 + (-t119 * t162 - t128 * t142 + (-t192 * t217 + t195 * t215) * t188 + (-t128 * t180 - t176 * t192) * t126) * t227 + (-t176 * t215 - t119 * t180 - t128 * t167 + t170 * t192 + (t128 * t162 - t192 * t220) * t126) * t228) * t134) * t129, (-t123 * t233 - t133 * t166) * t214 + (t123 * t236 + t145 * t133 + (-t166 * t122 - t123 * t144 + (-t121 * t162 - t127 * t142 + t168 + (-t127 * t180 - t164) * t126) * t227 + (-t121 * t180 - t127 * t167 - t143 + (t127 * t162 - t181) * t126) * t228) * t134) * t129, 0, 0, 0; 0, (-t147 * t158 + t159 * t230) * t213 + ((t171 * t189 + t203 * t186) * t147 + t159 * t237 + (-t158 * t141 - (-t171 * t186 + t203 * t189) * t150 - t159 * t140) * t148) * t137, -(-t148 * t229 + t231) * t204 * t213 + (t204 * t189 * t237 - t144 * t231 + (t144 * t229 - (t140 * t189 + t141 * t186) * t204) * t148) * t137, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:46
	% DurationCPUTime: 1.47s
	% Computational Cost: add. (4171->104), mult. (12384->230), div. (514->12), fcn. (16056->13), ass. (0->103)
	t208 = sin(qJ(2));
	t210 = cos(qJ(2));
	t258 = sin(pkin(10));
	t260 = cos(pkin(6));
	t230 = t260 * t258;
	t259 = cos(pkin(10));
	t200 = -t208 * t230 + t259 * t210;
	t204 = sin(pkin(11));
	t206 = cos(pkin(11));
	t207 = sin(qJ(3));
	t209 = cos(qJ(3));
	t231 = t260 * t259;
	t218 = -t208 * t231 - t258 * t210;
	t205 = sin(pkin(6));
	t237 = t205 * t259;
	t220 = t207 * t237 + t209 * t218;
	t229 = t258 * t208 - t210 * t231;
	t173 = -t220 * t204 - t229 * t206;
	t168 = t173 ^ 2;
	t244 = t208 * t209;
	t221 = t205 * t244 + t260 * t207;
	t245 = t206 * t210;
	t186 = t221 * t204 + t205 * t245;
	t183 = 0.1e1 / t186 ^ 2;
	t159 = t168 * t183 + 0.1e1;
	t187 = t207 * t218 - t209 * t237;
	t195 = t229 * qJD(2);
	t243 = qJD(2) * t206;
	t161 = (t187 * qJD(3) - t195 * t209) * t204 + t218 * t243;
	t235 = t260 * t209;
	t241 = qJD(3) * t207;
	t242 = qJD(2) * t210;
	t177 = qJD(3) * t204 * t235 + ((-t208 * t241 + t209 * t242) * t204 - t208 * t243) * t205;
	t182 = 0.1e1 / t186;
	t253 = t177 * t182 * t183;
	t254 = t173 * t183;
	t261 = -0.2e1 * (t161 * t254 - t168 * t253) / t159 ^ 2;
	t160 = atan2(-t173, t186);
	t152 = sin(t160);
	t153 = cos(t160);
	t227 = t152 * t173 - t153 * t186;
	t148 = 0.1e1 / t227;
	t217 = -t259 * t208 - t210 * t230;
	t236 = t205 * t258;
	t219 = -t200 * t209 - t207 * t236;
	t176 = -t204 * t217 - t206 * t219;
	t170 = 0.1e1 / t176;
	t149 = 0.1e1 / t227 ^ 2;
	t171 = 0.1e1 / t176 ^ 2;
	t156 = 0.1e1 / t159;
	t140 = (-t161 * t182 + t177 * t254) * t156;
	t228 = -t152 * t186 - t153 * t173;
	t137 = t228 * t140 - t152 * t161 + t153 * t177;
	t257 = t137 * t148 * t149;
	t175 = -t204 * t219 + t206 * t217;
	t256 = t149 * t175;
	t188 = -t200 * t207 + t209 * t236;
	t196 = t217 * qJD(2);
	t167 = t188 * qJD(3) + t196 * t209;
	t197 = t200 * qJD(2);
	t163 = t167 * t206 + t197 * t204;
	t172 = t170 * t171;
	t255 = t163 * t172;
	t249 = t217 * t209;
	t180 = t200 * t204 + t206 * t249;
	t252 = t180 * t188;
	t185 = t188 ^ 2;
	t251 = t185 * t206;
	t250 = t217 * t207;
	t248 = t204 * t209;
	t247 = t204 * t210;
	t246 = t205 * t207;
	t169 = t175 ^ 2;
	t145 = t149 * t169 + 0.1e1;
	t162 = t167 * t204 - t197 * t206;
	t240 = 0.2e1 * (t162 * t256 + t169 * t257) / t145 ^ 2;
	t158 = t171 * t185 + 0.1e1;
	t166 = t219 * qJD(3) - t196 * t207;
	t238 = t166 * t171 * t188;
	t239 = 0.2e1 * (-t185 * t255 + t238) / t158 ^ 2;
	t233 = -0.2e1 * t175 * t257;
	t232 = -0.2e1 * t173 * t253;
	t222 = t204 * t229;
	t178 = t206 * t218 - t209 * t222;
	t192 = (-t206 * t208 + t209 * t247) * t205;
	t225 = -t178 * t182 + t192 * t254;
	t201 = -t208 * t246 + t235;
	t224 = -t182 * t187 + t201 * t254;
	t223 = -t197 * t209 - t217 * t241;
	t190 = -t221 * qJD(3) - t242 * t246;
	t181 = (-t241 * t247 + (-t204 * t244 - t245) * qJD(2)) * t205;
	t179 = -t200 * t206 + t217 * t248;
	t165 = t220 * qJD(3) + t195 * t207;
	t164 = t218 * qJD(2) * t248 + t195 * t206 + t222 * t241;
	t154 = 0.1e1 / t158;
	t143 = 0.1e1 / t145;
	t142 = t224 * t204 * t156;
	t141 = t225 * t156;
	t139 = (-t152 * t187 + t153 * t201) * t204 + t228 * t142;
	t138 = t228 * t141 - t152 * t178 + t153 * t192;
	t136 = (t224 * t261 + (t201 * t232 - t165 * t182 + (t161 * t201 + t173 * t190 + t177 * t187) * t183) * t156) * t204;
	t135 = t225 * t261 + (t192 * t232 - t164 * t182 + (t161 * t192 + t173 * t181 + t177 * t178) * t183) * t156;
	t1 = [0, t135, t136, 0, 0, 0; 0, (t138 * t256 + t148 * t179) * t240 + (-(-t196 * t206 + t223 * t204) * t148 + t138 * t233 + (-t179 * t137 - t138 * t162 + (-(-t135 * t173 - t141 * t161 + t181 + (-t141 * t186 - t178) * t140) * t153 - (-t135 * t186 - t141 * t177 - t164 + (t141 * t173 - t192) * t140) * t152) * t175) * t149) * t143, (t148 * t188 * t204 + t139 * t256) * t240 + (-(t228 * t136 + (t227 * t140 - t152 * t177 - t153 * t161) * t142) * t256 + (-t149 * t162 + t233) * t139 + (-t166 * t148 + (-t188 * t137 - (-t152 * t165 + t153 * t190 + (-t152 * t201 - t153 * t187) * t140) * t175) * t149) * t204) * t143, 0, 0, 0; 0, (t170 * t250 + t171 * t252) * t239 + (0.2e1 * t252 * t255 + (-qJD(3) * t249 + t197 * t207) * t170 + (t163 * t250 - (t196 * t204 + t223 * t206) * t188 - t180 * t166) * t171) * t154, (-t170 * t219 + t171 * t251) * t239 + (-0.2e1 * t206 * t238 - t167 * t170 + (-t171 * t219 + 0.2e1 * t172 * t251) * t163) * t154, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:46
	% DurationCPUTime: 1.57s
	% Computational Cost: add. (3755->120), mult. (11172->251), div. (553->12), fcn. (14363->15), ass. (0->113)
	t258 = sin(pkin(10));
	t261 = cos(pkin(10));
	t265 = sin(qJ(2));
	t262 = cos(pkin(6));
	t268 = cos(qJ(2));
	t295 = t262 * t268;
	t248 = -t258 * t265 + t261 * t295;
	t241 = t248 * qJD(2);
	t296 = t262 * t265;
	t249 = t258 * t268 + t261 * t296;
	t264 = sin(qJ(3));
	t259 = sin(pkin(6));
	t299 = t259 * t264;
	t286 = t261 * t299;
	t267 = cos(qJ(3));
	t292 = qJD(3) * t267;
	t219 = -qJD(3) * t286 + t241 * t264 + t249 * t292;
	t298 = t259 * t267;
	t235 = t249 * t264 + t261 * t298;
	t233 = t235 ^ 2;
	t277 = -t262 * t267 + t265 * t299;
	t246 = 0.1e1 / t277 ^ 2;
	t229 = t233 * t246 + 0.1e1;
	t227 = 0.1e1 / t229;
	t253 = -t262 * t264 - t265 * t298;
	t293 = qJD(2) * t268;
	t285 = t259 * t293;
	t239 = qJD(3) * t253 - t264 * t285;
	t245 = 0.1e1 / t277;
	t304 = t235 * t246;
	t190 = (-t219 * t245 - t239 * t304) * t227;
	t230 = atan2(t235, -t277);
	t225 = sin(t230);
	t226 = cos(t230);
	t282 = t225 * t277 + t226 * t235;
	t185 = t282 * t190 + t225 * t219 + t226 * t239;
	t203 = t225 * t235 - t226 * t277;
	t200 = 0.1e1 / t203;
	t201 = 0.1e1 / t203 ^ 2;
	t315 = t185 * t200 * t201;
	t287 = t258 * t296;
	t251 = t261 * t268 - t287;
	t237 = -t251 * t264 + t258 * t298;
	t314 = 0.2e1 * t237 * t315;
	t302 = t239 * t245 * t246;
	t313 = (t219 * t304 + t233 * t302) / t229 ^ 2;
	t297 = t259 * t268;
	t288 = t235 * t297;
	t275 = -t245 * t248 + t246 * t288;
	t312 = t264 * t275;
	t257 = sin(pkin(11));
	t260 = cos(pkin(11));
	t263 = sin(qJ(6));
	t266 = cos(qJ(6));
	t280 = t257 * t266 - t260 * t263;
	t216 = t280 * t237;
	t238 = t251 * t267 + t258 * t299;
	t250 = t258 * t295 + t261 * t265;
	t223 = t238 * t257 - t250 * t260;
	t224 = t238 * t260 + t250 * t257;
	t209 = t223 * t263 + t224 * t266;
	t205 = 0.1e1 / t209;
	t206 = 0.1e1 / t209 ^ 2;
	t243 = t250 * qJD(2);
	t222 = qJD(3) * t237 - t243 * t267;
	t244 = -qJD(2) * t287 + t261 * t293;
	t212 = t222 * t257 - t244 * t260;
	t213 = t222 * t260 + t244 * t257;
	t188 = t209 * qJD(6) - t212 * t266 + t213 * t263;
	t283 = t223 * t266 - t224 * t263;
	t204 = t283 ^ 2;
	t193 = t204 * t206 + 0.1e1;
	t308 = t206 * t283;
	t189 = t283 * qJD(6) + t212 * t263 + t213 * t266;
	t310 = t189 * t205 * t206;
	t311 = (-t188 * t308 - t204 * t310) / t193 ^ 2;
	t309 = t201 * t237;
	t221 = -qJD(3) * t238 + t243 * t264;
	t307 = t221 * t201;
	t306 = t225 * t237;
	t305 = t226 * t237;
	t303 = t235 * t253;
	t301 = t250 * t264;
	t300 = t250 * t267;
	t294 = qJD(2) * t265;
	t234 = t237 ^ 2;
	t199 = t234 * t201 + 0.1e1;
	t291 = 0.2e1 * (-t234 * t315 + t237 * t307) / t199 ^ 2;
	t290 = 0.2e1 * t311;
	t284 = -0.2e1 * t283 * t310;
	t231 = -t251 * t260 - t257 * t300;
	t232 = t251 * t257 - t260 * t300;
	t281 = t231 * t266 - t232 * t263;
	t211 = t231 * t263 + t232 * t266;
	t279 = t257 * t263 + t260 * t266;
	t236 = t249 * t267 - t286;
	t278 = t236 * t245 + t246 * t303;
	t276 = qJD(3) * t301 - t244 * t267;
	t217 = t279 * t237;
	t242 = t249 * qJD(2);
	t240 = qJD(3) * t277 - t267 * t285;
	t220 = -t235 * qJD(3) + t241 * t267;
	t215 = -t243 * t257 + t260 * t276;
	t214 = t243 * t260 + t257 * t276;
	t196 = 0.1e1 / t199;
	t195 = t227 * t312;
	t194 = t278 * t227;
	t191 = 0.1e1 / t193;
	t187 = (t225 * t248 - t226 * t297) * t264 + t282 * t195;
	t186 = -t282 * t194 + t225 * t236 + t226 * t253;
	t183 = 0.2e1 * t278 * t313 + (-0.2e1 * t302 * t303 - t220 * t245 + (-t219 * t253 - t235 * t240 - t236 * t239) * t246) * t227;
	t181 = -0.2e1 * t312 * t313 + (t275 * t292 + (0.2e1 * t288 * t302 + t242 * t245 + (-t239 * t248 + (t219 * t268 - t235 * t294) * t259) * t246) * t264) * t227;
	t1 = [0, t181, t183, 0, 0, 0; 0, (t187 * t309 - t200 * t301) * t291 + ((t244 * t264 + t250 * t292) * t200 + (-t307 + t314) * t187 + (-t301 * t185 - (t181 * t235 + t195 * t219 + (t264 * t294 - t268 * t292) * t259 + (t195 * t277 + t248 * t264) * t190) * t305 - (t248 * t292 + t181 * t277 - t195 * t239 - t242 * t264 + (-t195 * t235 + t264 * t297) * t190) * t306) * t201) * t196, (t186 * t309 + t200 * t238) * t291 + (t186 * t314 - t222 * t200 + (t238 * t185 - t186 * t221 - (t183 * t235 - t194 * t219 + t240 + (-t194 * t277 + t236) * t190) * t305 - (t183 * t277 + t194 * t239 + t220 + (t194 * t235 - t253) * t190) * t306) * t201) * t196, 0, 0, 0; 0, (t205 * t281 - t211 * t308) * t290 + ((t211 * qJD(6) - t214 * t266 + t215 * t263) * t205 + t211 * t284 + (t281 * t189 + (t281 * qJD(6) + t214 * t263 + t215 * t266) * t283 - t211 * t188) * t206) * t191, (t205 * t216 - t217 * t308) * t290 + ((qJD(6) * t217 - t280 * t221) * t205 + t217 * t284 + (t216 * t189 + (qJD(6) * t216 + t279 * t221) * t283 - t217 * t188) * t206) * t191, 0, 0, -0.2e1 * t311 - 0.2e1 * (t188 * t206 * t191 - (-t191 * t310 - t206 * t311) * t283) * t283;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end