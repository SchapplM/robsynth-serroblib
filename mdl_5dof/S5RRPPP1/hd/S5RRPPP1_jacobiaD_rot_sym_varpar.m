% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPPP1
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
%   Wie in S5RRPPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRPPP1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPP1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:25:39
	% EndTime: 2019-12-31 19:25:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:25:38
	% EndTime: 2019-12-31 19:25:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:25:39
	% EndTime: 2019-12-31 19:25:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:25:39
	% EndTime: 2019-12-31 19:25:40
	% DurationCPUTime: 0.74s
	% Computational Cost: add. (1211->84), mult. (3693->201), div. (600->14), fcn. (4788->11), ass. (0->93)
	t139 = sin(pkin(8));
	t141 = cos(pkin(8));
	t143 = sin(qJ(2));
	t142 = cos(pkin(5));
	t146 = cos(qJ(1));
	t180 = t146 * t142;
	t140 = sin(pkin(5));
	t144 = sin(qJ(1));
	t188 = t140 * t144;
	t156 = t143 * t180 - t188;
	t145 = cos(qJ(2));
	t181 = t145 * t146;
	t115 = -t156 * t139 + t141 * t181;
	t114 = t139 * t181 + t156 * t141;
	t184 = t144 * t142;
	t187 = t140 * t146;
	t130 = t143 * t187 + t184;
	t182 = t145 * t140;
	t162 = qJD(2) * t182;
	t117 = t130 * qJD(1) + t144 * t162;
	t166 = t143 * t188;
	t128 = t166 - t180;
	t134 = 0.1e1 / t140;
	t136 = 0.1e1 / t145;
	t137 = 0.1e1 / t145 ^ 2;
	t176 = qJD(2) * t143;
	t163 = t137 * t176;
	t202 = (t117 * t136 + t128 * t163) * t134;
	t120 = atan2(-t128, -t182);
	t118 = sin(t120);
	t119 = cos(t120);
	t126 = t128 ^ 2;
	t135 = 0.1e1 / t140 ^ 2;
	t123 = t126 * t135 * t137 + 0.1e1;
	t121 = 0.1e1 / t123;
	t190 = t128 * t134;
	t168 = t136 * t190;
	t199 = (t119 * t168 - t118) * t121 + t118;
	t192 = t118 * t128;
	t102 = -t119 * t182 - t192;
	t98 = 0.1e1 / t102;
	t109 = 0.1e1 / t115;
	t99 = 0.1e1 / t102 ^ 2;
	t110 = 0.1e1 / t115 ^ 2;
	t177 = qJD(1) * t146;
	t116 = qJD(1) * t166 - t142 * t177 - t146 * t162;
	t127 = t130 ^ 2;
	t194 = t130 * t99;
	t93 = t121 * t202;
	t89 = (-t128 * t93 + t140 * t176) * t119 + (t93 * t182 - t117) * t118;
	t196 = t98 * t99 * t89;
	t96 = t127 * t99 + 0.1e1;
	t198 = (-t116 * t194 - t127 * t196) / t96 ^ 2;
	t94 = 0.1e1 / t96;
	t197 = t94 * t99;
	t138 = t136 * t137;
	t195 = 0.1e1 / t123 ^ 2 * (t117 * t128 * t137 + t126 * t138 * t176) * t135;
	t193 = t110 * t114;
	t191 = t119 * t128;
	t189 = t137 * t143;
	t186 = t142 * t143;
	t185 = t142 * t145;
	t183 = t144 * t145;
	t159 = t189 * t190 + t144;
	t101 = t159 * t121;
	t179 = t101 - t144;
	t178 = qJD(1) * t144;
	t175 = qJD(2) * t146;
	t108 = t114 ^ 2;
	t105 = t108 * t110 + 0.1e1;
	t155 = t143 * t184 + t187;
	t112 = -t139 * t183 - t155 * t141;
	t158 = -t139 * t143 + t141 * t185;
	t124 = t158 * t146;
	t106 = t112 * qJD(1) + qJD(2) * t124;
	t113 = t155 * t139 - t141 * t183;
	t157 = t139 * t185 + t141 * t143;
	t125 = t157 * t146;
	t107 = t113 * qJD(1) - qJD(2) * t125;
	t111 = t109 * t110;
	t174 = 0.2e1 / t105 ^ 2 * (-t107 * t108 * t111 + t106 * t193);
	t173 = -0.2e1 * t195;
	t172 = t98 * t198;
	t171 = 0.2e1 * t111 * t114;
	t170 = t116 * t197;
	t169 = t94 * t194;
	t161 = t118 * t169;
	t153 = 0.2e1 * t94 * t196 + 0.2e1 * t99 * t198;
	t103 = 0.1e1 / t105;
	t92 = t199 * t130;
	t90 = -t101 * t191 + (t179 * t145 * t118 + t119 * t143) * t140;
	t88 = t159 * t173 + (t177 + (t117 * t189 + (0.2e1 * t138 * t143 ^ 2 + t136) * t128 * qJD(2)) * t134) * t121;
	t1 = [(t130 * t136 * t173 + (-t116 * t136 + t130 * t163) * t121) * t134, t88, 0, 0, 0; 0.2e1 * t128 * t172 + (-t117 * t98 + (-t116 * t92 + t128 * t89) * t99) * t94 + (-t92 * t153 - (t121 * t93 * t168 + t173) * t161 - (0.2e1 * t168 * t195 - t93 + (t93 - t202) * t121) * t119 * t169 - t199 * t170) * t130, t90 * t170 + (-(-t88 * t191 + (-t117 * t119 + t93 * t192) * t101) * t197 + t90 * t153) * t130 + ((-t98 * t94 * t175 - (-t179 * qJD(2) - t93) * t161) * t143 + (-0.2e1 * t146 * t172 + (-t98 * t178 + (-t146 * t89 + (-(t88 - t177) * t118 - (t179 * t93 + qJD(2)) * t119) * t130) * t99) * t94) * t145) * t140, 0, 0, 0; (-t109 * t112 + t113 * t193) * t174 + (-t113 * t106 * t110 + (-t112 * t110 + t113 * t171) * t107 + (-t158 * t109 - t157 * t193) * t144 * qJD(2) + (-t114 * t109 + t115 * t193) * qJD(1)) * t103, (-t109 * t124 - t125 * t193) * t174 + (((-t139 * t145 - t141 * t186) * t175 - t158 * t178) * t109 - t125 * t107 * t171 + (-t124 * t107 - ((t139 * t186 - t141 * t145) * t175 + t157 * t178) * t114 + t125 * t106) * t110) * t103, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:25:39
	% EndTime: 2019-12-31 19:25:40
	% DurationCPUTime: 0.84s
	% Computational Cost: add. (2054->91), mult. (6960->199), div. (427->12), fcn. (8805->11), ass. (0->96)
	t162 = sin(pkin(8));
	t165 = cos(pkin(5));
	t168 = cos(qJ(2));
	t164 = cos(pkin(8));
	t166 = sin(qJ(2));
	t209 = t164 * t166;
	t154 = t162 * t168 + t165 * t209;
	t146 = t154 * qJD(2);
	t208 = t165 * t168;
	t194 = t164 * t208;
	t181 = -t162 * t166 + t194;
	t148 = 0.1e1 / t181 ^ 2;
	t214 = t146 * t148;
	t167 = sin(qJ(1));
	t205 = t167 * t165;
	t163 = sin(pkin(5));
	t169 = cos(qJ(1));
	t210 = t163 * t169;
	t178 = t166 * t205 + t210;
	t204 = t167 * t168;
	t136 = t162 * t204 + t178 * t164;
	t130 = atan2(-t136, -t181);
	t123 = sin(t130);
	t124 = cos(t130);
	t133 = t136 ^ 2;
	t129 = t133 * t148 + 0.1e1;
	t125 = 0.1e1 / t129;
	t147 = 0.1e1 / t181;
	t218 = t136 * t147;
	t225 = (-t124 * t218 + t123) * t125 - t123;
	t203 = t168 * t169;
	t206 = t166 * t169;
	t211 = t163 * t167;
	t224 = (t165 * t206 - t211) * t162 - t164 * t203;
	t119 = -t123 * t136 - t124 * t181;
	t116 = 0.1e1 / t119;
	t156 = t163 * t206 + t205;
	t150 = 0.1e1 / t156;
	t117 = 0.1e1 / t119 ^ 2;
	t151 = 0.1e1 / t156 ^ 2;
	t177 = t154 * t169;
	t186 = qJD(2) * t194;
	t202 = qJD(1) * t167;
	t190 = t164 * t202;
	t207 = t166 * t167;
	t195 = t162 * t207;
	t122 = qJD(1) * t177 - qJD(2) * t195 - t163 * t190 + t167 * t186;
	t183 = t122 * t147 + t136 * t214;
	t109 = t183 * t125;
	t184 = t123 * t181 - t124 * t136;
	t106 = t184 * t109 - t123 * t122 + t124 * t146;
	t223 = t106 * t116 * t117;
	t213 = t147 * t214;
	t222 = (t122 * t136 * t148 + t133 * t213) / t129 ^ 2;
	t139 = -t164 * t211 + t177;
	t221 = t117 * t139;
	t220 = t123 * t139;
	t219 = t124 * t139;
	t217 = t136 * t154;
	t216 = t224 * t151;
	t155 = -t163 * t207 + t165 * t169;
	t189 = qJD(2) * t163 * t168;
	t141 = t155 * qJD(1) + t169 * t189;
	t215 = t141 * t150 * t151;
	t212 = t162 * t165;
	t201 = qJD(2) * t169;
	t134 = t139 ^ 2;
	t114 = t117 * t134 + 0.1e1;
	t188 = t166 * t201;
	t120 = t136 * qJD(1) + t162 * t188 - t169 * t186;
	t200 = 0.2e1 * (-t120 * t221 - t134 * t223) / t114 ^ 2;
	t199 = 0.2e1 * t223;
	t192 = t164 * t204;
	t138 = t178 * t162 - t192;
	t180 = t162 * t208 + t209;
	t144 = t180 * t169;
	t121 = t138 * qJD(1) - qJD(2) * t144;
	t135 = t224 ^ 2;
	t131 = t135 * t151 + 0.1e1;
	t198 = 0.2e1 * (-t121 * t216 - t135 * t215) / t131 ^ 2;
	t197 = -0.2e1 * t222;
	t196 = t147 * t222;
	t187 = -0.2e1 * t224 * t215;
	t142 = -t165 * t192 + t195;
	t182 = -t142 * t147 + t148 * t217;
	t176 = qJD(1) * t181;
	t145 = t181 * qJD(2);
	t143 = t181 * t169;
	t132 = t167 * t146 - t169 * t176;
	t127 = 0.1e1 / t131;
	t112 = 0.1e1 / t114;
	t110 = t182 * t125;
	t108 = t225 * t139;
	t107 = t184 * t110 + t123 * t142 + t124 * t154;
	t105 = t182 * t197 + (0.2e1 * t213 * t217 - t132 * t147 + (t122 * t154 + t136 * t145 - t142 * t146) * t148) * t125;
	t1 = [-0.2e1 * t139 * t196 + (-t120 * t147 + t139 * t214) * t125, t105, 0, 0, 0; t136 * t116 * t200 + (-t122 * t116 + (t106 * t136 + t108 * t120) * t117) * t112 + (t108 * t199 * t112 + (t108 * t200 + (-(t109 * t125 * t218 + t197) * t220 - (0.2e1 * t136 * t196 - t109 + (t109 - t183) * t125) * t219 + t225 * t120) * t112) * t117) * t139, (t107 * t221 - t116 * t143) * t200 + ((-t154 * t201 - t167 * t176) * t116 + t107 * t139 * t199 + (-t143 * t106 + t107 * t120 - (-t105 * t136 - t110 * t122 + t145 + (t110 * t181 + t142) * t109) * t219 - (t105 * t181 - t110 * t146 + t132 + (t110 * t136 - t154) * t109) * t220) * t117) * t112, 0, 0, 0; (-t138 * t150 - t155 * t216) * t198 + (t155 * t187 + (-t138 * t141 + (-t156 * qJD(1) - t167 * t189) * t224 - t155 * t121) * t151 + (t180 * t167 * qJD(2) + t224 * qJD(1)) * t150) * t127, (-t163 * t203 * t216 + t144 * t150) * t198 + ((t201 * t212 + t190) * t150 * t166 + ((-t164 * t201 + t202 * t212) * t150 + t187 * t210) * t168 + (t144 * t141 + (-t121 * t203 - (t168 * t202 + t188) * t224) * t163) * t151) * t127, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:25:39
	% EndTime: 2019-12-31 19:25:40
	% DurationCPUTime: 0.83s
	% Computational Cost: add. (2054->88), mult. (6960->192), div. (427->12), fcn. (8805->11), ass. (0->101)
	t179 = sin(qJ(1));
	t177 = cos(pkin(5));
	t178 = sin(qJ(2));
	t244 = sin(pkin(8));
	t209 = t244 * t178;
	t204 = t177 * t209;
	t176 = cos(pkin(8));
	t180 = cos(qJ(2));
	t223 = t180 * t176;
	t192 = t204 - t223;
	t181 = cos(qJ(1));
	t175 = sin(pkin(5));
	t210 = t175 * t244;
	t202 = t181 * t210;
	t147 = -t192 * t179 - t202;
	t251 = 0.2e1 * t147;
	t143 = t147 ^ 2;
	t208 = t244 * t180;
	t226 = t178 * t176;
	t191 = -t177 * t208 - t226;
	t158 = 0.1e1 / t191 ^ 2;
	t139 = t143 * t158 + 0.1e1;
	t135 = 0.1e1 / t139;
	t152 = t191 * t179;
	t222 = t180 * t181;
	t190 = t176 * t222 + t179 * t210;
	t197 = qJD(1) * t204;
	t132 = qJD(1) * t190 + qJD(2) * t152 - t181 * t197;
	t157 = 0.1e1 / t191;
	t156 = t192 * qJD(2);
	t230 = t156 * t158;
	t196 = t132 * t157 - t147 * t230;
	t119 = t196 * t135;
	t140 = atan2(-t147, -t191);
	t133 = sin(t140);
	t134 = cos(t140);
	t198 = t133 * t191 - t134 * t147;
	t116 = t119 * t198 - t132 * t133 - t134 * t156;
	t129 = -t133 * t147 - t134 * t191;
	t127 = 0.1e1 / t129 ^ 2;
	t250 = t116 * t127;
	t228 = t175 * t179;
	t165 = t177 * t181 - t178 * t228;
	t211 = qJD(2) * t175 * t180;
	t151 = qJD(1) * t165 + t181 * t211;
	t224 = t179 * t177;
	t225 = t178 * t181;
	t166 = t175 * t225 + t224;
	t160 = 0.1e1 / t166;
	t247 = t176 * (t177 * t225 - t228) + t181 * t208;
	t144 = t247 ^ 2;
	t161 = 0.1e1 / t166 ^ 2;
	t141 = t144 * t161 + 0.1e1;
	t137 = 0.1e1 / t141;
	t235 = t137 * t161;
	t227 = t175 * t181;
	t146 = t179 * t208 + (t178 * t224 + t227) * t176;
	t193 = t177 * t223 - t209;
	t153 = t193 * t181;
	t130 = qJD(1) * t146 - qJD(2) * t153;
	t231 = t151 * t160 * t161;
	t232 = t247 * t161;
	t240 = (-t130 * t232 - t144 * t231) / t141 ^ 2;
	t249 = -t151 * t235 - 0.2e1 * t160 * t240;
	t234 = t147 * t157;
	t248 = t135 * (-t134 * t234 + t133) - t133;
	t246 = -0.2e1 * t247 * (t137 * t231 + t161 * t240);
	t126 = 0.1e1 / t129;
	t150 = -t181 * t204 + t190;
	t145 = t150 ^ 2;
	t124 = t127 * t145 + 0.1e1;
	t212 = qJD(1) * t179 * t180;
	t220 = qJD(2) * t181;
	t131 = -qJD(1) * t202 + t176 * t212 - t179 * t197 - t191 * t220;
	t237 = t131 * t127;
	t242 = t126 * t250;
	t243 = (-t145 * t242 - t150 * t237) / t124 ^ 2;
	t233 = t147 * t192;
	t195 = t152 * t157 - t158 * t233;
	t120 = t195 * t135;
	t117 = t120 * t198 - t133 * t152 - t134 * t192;
	t122 = 0.1e1 / t124;
	t241 = t117 * t122;
	t229 = t157 * t230;
	t239 = (t132 * t147 * t158 - t143 * t229) / t139 ^ 2;
	t238 = t127 * t150;
	t236 = t137 * t160;
	t219 = 0.2e1 * t242;
	t218 = -0.2e1 * t239;
	t216 = t126 * t243;
	t215 = t157 * t239;
	t214 = t137 * t232;
	t207 = 0.2e1 * t127 * t243;
	t189 = t193 * t179;
	t188 = qJD(1) * t191;
	t155 = t191 * qJD(2);
	t154 = t191 * t181;
	t142 = -t156 * t179 - t181 * t188;
	t118 = t248 * t150;
	t115 = t195 * t218 + (0.2e1 * t229 * t233 - t142 * t157 + (-t132 * t192 + t147 * t155 - t152 * t156) * t158) * t135;
	t1 = [-0.2e1 * t150 * t215 + (-t131 * t157 - t150 * t230) * t135, t115, 0, 0, 0; t216 * t251 + (-t132 * t126 + (t116 * t147 + t118 * t131) * t127) * t122 + (t118 * t207 + (t118 * t219 + t248 * t237 + (-(t119 * t135 * t234 + t218) * t133 - (t215 * t251 - t119 + (t119 - t196) * t135) * t134) * t238) * t122) * t150, -0.2e1 * t154 * t216 + t237 * t241 + (t117 * t207 + t219 * t241) * t150 + ((-t179 * t188 + t192 * t220) * t126 - t154 * t250 - ((-t115 * t147 - t120 * t132 + t155 + (t120 * t191 - t152) * t119) * t134 + (t115 * t191 + t120 * t156 + t142 + (t120 * t147 + t192) * t119) * t133) * t238) * t122, 0, 0, 0; (t247 * qJD(1) + qJD(2) * t189) * t236 + (-qJD(1) * t166 - t179 * t211) * t214 + t249 * t146 + (-t130 * t235 + t246) * t165, ((t177 * t226 + t208) * t220 + qJD(1) * t189) * t236 + (-qJD(2) * t178 * t247 - t180 * t130) * t227 * t235 - t249 * t153 + (-t212 * t214 + t222 * t246) * t175, 0, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end