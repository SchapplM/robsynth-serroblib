% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR6
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
%   Wie in S6RRRPRR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:02
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:22
	% DurationCPUTime: 1.01s
	% Computational Cost: add. (1002->94), mult. (2519->211), div. (480->12), fcn. (2968->9), ass. (0->92)
	t126 = sin(qJ(1));
	t119 = t126 ^ 2;
	t125 = sin(qJ(2));
	t118 = t125 ^ 2;
	t128 = cos(qJ(2));
	t121 = 0.1e1 / t128 ^ 2;
	t173 = t118 * t121;
	t113 = t119 * t173 + 0.1e1;
	t117 = t125 * t118;
	t120 = 0.1e1 / t128;
	t172 = t120 * t125;
	t137 = qJD(2) * (t117 * t120 * t121 + t172);
	t129 = cos(qJ(1));
	t163 = qJD(1) * t129;
	t151 = t126 * t163;
	t181 = 0.1e1 / t113 ^ 2 * (t119 * t137 + t151 * t173);
	t193 = -0.2e1 * t181;
	t111 = 0.1e1 / t113;
	t146 = 0.1e1 + t173;
	t190 = t126 * t146;
	t98 = t111 * t190;
	t192 = t126 * t98 - 0.1e1;
	t124 = sin(qJ(3));
	t127 = cos(qJ(3));
	t165 = t129 * t127;
	t107 = t126 * t124 + t128 * t165;
	t102 = 0.1e1 / t107 ^ 2;
	t166 = t129 * t124;
	t168 = t126 * t127;
	t106 = t128 * t166 - t168;
	t175 = t106 * t127;
	t101 = 0.1e1 / t107;
	t177 = t101 * t124;
	t139 = t175 * t102 - t177;
	t100 = t106 ^ 2;
	t99 = t100 * t102 + 0.1e1;
	t96 = 0.1e1 / t99;
	t191 = t139 * t96;
	t169 = t126 * t125;
	t110 = atan2(-t169, -t128);
	t109 = cos(t110);
	t108 = sin(t110);
	t154 = t108 * t169;
	t94 = -t109 * t128 - t154;
	t91 = 0.1e1 / t94;
	t92 = 0.1e1 / t94 ^ 2;
	t189 = t111 - 0.1e1;
	t123 = t129 ^ 2;
	t161 = qJD(2) * t128;
	t155 = t92 * t161;
	t150 = t125 * t163;
	t162 = qJD(2) * t126;
	t174 = t109 * t125;
	t149 = t121 * t162;
	t85 = (-(-t126 * t161 - t150) * t120 + t118 * t149) * t111;
	t80 = (-t126 * t85 + qJD(2)) * t174 + (-t150 + (t85 - t162) * t128) * t108;
	t187 = t80 * t91 * t92;
	t90 = t123 * t118 * t92 + 0.1e1;
	t188 = (t123 * t125 * t155 + (-t123 * t187 - t151 * t92) * t118) / t90 ^ 2;
	t176 = t102 * t106;
	t143 = -qJD(1) * t128 + qJD(3);
	t144 = qJD(3) * t128 - qJD(1);
	t160 = qJD(2) * t129;
	t148 = t125 * t160;
	t87 = -t144 * t166 + (t126 * t143 - t148) * t127;
	t183 = t101 * t102 * t87;
	t167 = t126 * t128;
	t138 = t124 * t167 + t165;
	t86 = qJD(1) * t138 - t107 * qJD(3) + t124 * t148;
	t186 = (-t100 * t183 - t176 * t86) / t99 ^ 2;
	t88 = 0.1e1 / t90;
	t185 = t88 * t92;
	t184 = t91 * t88;
	t179 = t129 * t92;
	t178 = qJD(2) * t98;
	t171 = t125 * t129;
	t164 = qJD(1) * t126;
	t159 = 0.2e1 * t187;
	t158 = -0.2e1 * t186;
	t157 = t91 * t188;
	t156 = t106 * t183;
	t153 = t111 * t118 * t120;
	t147 = 0.2e1 * t92 * t188;
	t145 = t120 * t193;
	t142 = t126 * t153;
	t141 = t146 * t129;
	t140 = t143 * t129;
	t105 = -t127 * t167 + t166;
	t84 = (t108 * t125 * t189 - t109 * t142) * t129;
	t83 = -t192 * t174 + (-t126 + t98) * t128 * t108;
	t81 = t190 * t193 + (qJD(1) * t141 + 0.2e1 * t126 * t137) * t111;
	t1 = [t145 * t171 + (qJD(2) * t141 - t164 * t172) * t111, t81, 0, 0, 0, 0; (-t161 * t184 + (0.2e1 * t157 + (qJD(1) * t84 + t80) * t185) * t125) * t126 + (t84 * t147 * t125 + (-t84 * t155 + (t84 * t159 + ((0.2e1 * t125 * t181 - t142 * t85 - t161 * t189) * t108 + (t118 * t126 * t145 + t125 * t85 + (t117 * t149 - (t85 - 0.2e1 * t162) * t125) * t111) * t109) * t179) * t125 + (-t91 + (-(t119 - t123) * t109 * t153 + t189 * t154) * t92) * t125 * qJD(1)) * t88) * t129, (-t164 * t184 + (-0.2e1 * t157 + (-qJD(2) * t83 - t80) * t185) * t129) * t128 + (t83 * t129 * t147 + (-t91 * t160 - ((-t126 * t81 - t163 * t98) * t109 + (t192 * t85 + t162 - t178) * t108) * t92 * t171 + (t129 * t159 + t164 * t92) * t83 - ((t81 - t163) * t108 + (t85 * t98 + qJD(2) + (-t85 - t178) * t126) * t109) * t128 * t179) * t88) * t125, 0, 0, 0, 0; 0.2e1 * (t101 * t138 + t105 * t176) * t186 + (0.2e1 * t105 * t156 - t144 * t101 * t168 + (t162 * t125 + t140) * t177 + (t138 * t87 + t105 * t86 - t140 * t175 - (qJD(2) * t125 * t127 + t124 * t144) * t106 * t126) * t102) * t96, t128 * t160 * t191 + (-t164 * t191 + (t139 * t158 + ((-qJD(3) * t101 - 0.2e1 * t156) * t127 + (-t127 * t86 + (-qJD(3) * t106 + t87) * t124) * t102) * t96) * t129) * t125, t158 + 0.2e1 * (-t86 * t102 * t96 + (-t102 * t186 - t183 * t96) * t106) * t106, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:22
	% DurationCPUTime: 1.03s
	% Computational Cost: add. (1350->93), mult. (2519->209), div. (480->12), fcn. (2968->9), ass. (0->94)
	t135 = sin(qJ(1));
	t129 = t135 ^ 2;
	t134 = sin(qJ(2));
	t128 = t134 ^ 2;
	t136 = cos(qJ(2));
	t131 = 0.1e1 / t136 ^ 2;
	t184 = t128 * t131;
	t121 = t129 * t184 + 0.1e1;
	t127 = t134 * t128;
	t130 = 0.1e1 / t136;
	t181 = t130 * t134;
	t146 = qJD(2) * (t127 * t130 * t131 + t181);
	t137 = cos(qJ(1));
	t173 = qJD(1) * t137;
	t182 = t128 * t135;
	t150 = t173 * t182;
	t190 = (t129 * t146 + t131 * t150) / t121 ^ 2;
	t200 = -0.2e1 * t190;
	t126 = qJ(3) + pkin(11);
	t125 = cos(t126);
	t175 = t136 * t137;
	t124 = sin(t126);
	t179 = t135 * t124;
	t114 = t125 * t175 + t179;
	t109 = 0.1e1 / t114 ^ 2;
	t178 = t135 * t125;
	t113 = t124 * t175 - t178;
	t187 = t113 * t125;
	t108 = 0.1e1 / t114;
	t189 = t108 * t124;
	t148 = t109 * t187 - t189;
	t107 = t113 ^ 2;
	t100 = t107 * t109 + 0.1e1;
	t98 = 0.1e1 / t100;
	t199 = t148 * t98;
	t157 = 0.1e1 + t184;
	t198 = t135 * t157;
	t177 = t135 * t134;
	t118 = atan2(-t177, -t136);
	t116 = cos(t118);
	t115 = sin(t118);
	t163 = t115 * t177;
	t104 = -t116 * t136 - t163;
	t101 = 0.1e1 / t104;
	t102 = 0.1e1 / t104 ^ 2;
	t119 = 0.1e1 / t121;
	t197 = t119 - 0.1e1;
	t133 = t137 ^ 2;
	t171 = qJD(2) * t136;
	t183 = t128 * t133;
	t172 = qJD(2) * t135;
	t159 = t131 * t172;
	t160 = t134 * t173;
	t94 = (-(-t135 * t171 - t160) * t130 + t128 * t159) * t119;
	t154 = t94 - t172;
	t155 = -t135 * t94 + qJD(2);
	t186 = t116 * t134;
	t88 = t155 * t186 + (t154 * t136 - t160) * t115;
	t193 = t101 * t102 * t88;
	t97 = t102 * t183 + 0.1e1;
	t196 = (-t183 * t193 + (t133 * t134 * t171 - t150) * t102) / t97 ^ 2;
	t188 = t109 * t113;
	t152 = -qJD(1) * t136 + qJD(3);
	t153 = qJD(3) * t136 - qJD(1);
	t170 = qJD(2) * t137;
	t158 = t134 * t170;
	t185 = t124 * t137;
	t93 = -t153 * t185 + (t152 * t135 - t158) * t125;
	t192 = t108 * t109 * t93;
	t176 = t135 * t136;
	t147 = t124 * t176 + t125 * t137;
	t92 = t147 * qJD(1) - qJD(3) * t114 + t124 * t158;
	t195 = (-t107 * t192 - t92 * t188) / t100 ^ 2;
	t95 = 0.1e1 / t97;
	t194 = t102 * t95;
	t180 = t134 * t137;
	t174 = qJD(1) * t135;
	t169 = 0.2e1 * t196;
	t168 = -0.2e1 * t195;
	t167 = 0.2e1 * t193;
	t166 = t101 * t196;
	t165 = t113 * t192;
	t164 = t95 * t171;
	t162 = t119 * t128 * t130;
	t156 = t130 * t200;
	t151 = t135 * t162;
	t149 = t157 * t137;
	t145 = t134 * t172 + t152 * t137;
	t112 = -t125 * t176 + t185;
	t106 = t119 * t198;
	t91 = (t197 * t134 * t115 - t116 * t151) * t137;
	t90 = -t115 * t176 + t186 + (t115 * t136 - t116 * t177) * t106;
	t89 = t198 * t200 + (qJD(1) * t149 + 0.2e1 * t135 * t146) * t119;
	t1 = [t156 * t180 + (qJD(2) * t149 - t174 * t181) * t119, t89, 0, 0, 0, 0; (-t101 * t164 + (0.2e1 * t166 + (qJD(1) * t91 + t88) * t194) * t134) * t135 + ((-t91 * t164 + (t91 * t169 + ((0.2e1 * t134 * t190 - t94 * t151 - t197 * t171) * t115 + (t156 * t182 + t134 * t94 + (t127 * t159 - (t94 - 0.2e1 * t172) * t134) * t119) * t116) * t95 * t137) * t134) * t102 + (t91 * t167 + (-t101 + ((-t129 + t133) * t116 * t162 + t197 * t163) * t102) * qJD(1)) * t134 * t95) * t137, (-t101 * t95 * t174 + (-0.2e1 * t166 + (-qJD(2) * t90 - t88) * t194) * t137) * t136 + (t90 * t137 * t102 * t169 + ((-qJD(2) * t101 + t90 * t167) * t137 + (t90 * t174 + (-(-t106 * t173 - t135 * t89) * t116 - ((t106 * t135 - 0.1e1) * t94 + (-t106 + t135) * qJD(2)) * t115) * t180) * t102) * t95 - ((t89 - t173) * t115 + (t154 * t106 + t155) * t116) * t175 * t194) * t134, 0, 0, 0, 0; 0.2e1 * (t108 * t147 + t112 * t188) * t195 + (0.2e1 * t112 * t165 - t153 * t108 * t178 + t145 * t189 + (-t153 * t113 * t179 + t112 * t92 - t145 * t187 + t147 * t93) * t109) * t98, t136 * t170 * t199 + (-t174 * t199 + (t148 * t168 + ((-qJD(3) * t108 - 0.2e1 * t165) * t125 + (-t125 * t92 + (-qJD(3) * t113 + t93) * t124) * t109) * t98) * t137) * t134, t168 + 0.2e1 * (-t109 * t92 * t98 + (-t109 * t195 - t98 * t192) * t113) * t113, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:22
	% DurationCPUTime: 1.09s
	% Computational Cost: add. (2009->96), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->95)
	t155 = sin(qJ(2));
	t149 = t155 ^ 2;
	t157 = cos(qJ(2));
	t152 = 0.1e1 / t157 ^ 2;
	t202 = t149 * t152;
	t156 = sin(qJ(1));
	t220 = 0.2e1 * t156;
	t219 = t155 * t202;
	t146 = qJ(3) + pkin(11) + qJ(5);
	t145 = cos(t146);
	t158 = cos(qJ(1));
	t194 = t157 * t158;
	t144 = sin(t146);
	t198 = t156 * t144;
	t134 = t145 * t194 + t198;
	t196 = t156 * t155;
	t139 = atan2(-t196, -t157);
	t137 = cos(t139);
	t136 = sin(t139);
	t183 = t136 * t196;
	t124 = -t137 * t157 - t183;
	t121 = 0.1e1 / t124;
	t128 = 0.1e1 / t134;
	t151 = 0.1e1 / t157;
	t122 = 0.1e1 / t124 ^ 2;
	t129 = 0.1e1 / t134 ^ 2;
	t218 = -0.2e1 * t155;
	t150 = t156 ^ 2;
	t142 = t150 * t202 + 0.1e1;
	t140 = 0.1e1 / t142;
	t217 = t140 - 0.1e1;
	t147 = qJD(3) + qJD(5);
	t195 = t156 * t157;
	t168 = t144 * t195 + t145 * t158;
	t189 = qJD(2) * t158;
	t179 = t155 * t189;
	t112 = t168 * qJD(1) - t134 * t147 + t144 * t179;
	t197 = t156 * t145;
	t133 = t144 * t194 - t197;
	t127 = t133 ^ 2;
	t117 = t127 * t129 + 0.1e1;
	t207 = t129 * t133;
	t173 = -qJD(1) * t157 + t147;
	t174 = t147 * t157 - qJD(1);
	t204 = t144 * t158;
	t113 = -t174 * t204 + (t173 * t156 - t179) * t145;
	t214 = t113 * t128 * t129;
	t216 = (-t112 * t207 - t127 * t214) / t117 ^ 2;
	t192 = qJD(1) * t158;
	t180 = t155 * t192;
	t190 = qJD(2) * t157;
	t191 = qJD(2) * t156;
	t114 = (-(-t156 * t190 - t180) * t151 + t191 * t202) * t140;
	t205 = t137 * t155;
	t108 = (-t114 * t156 + qJD(2)) * t205 + (-t180 + (t114 - t191) * t157) * t136;
	t215 = t108 * t121 * t122;
	t213 = t114 * t136;
	t212 = t114 * t155;
	t211 = t122 * t155;
	t200 = t151 * t155;
	t167 = qJD(2) * (t151 * t219 + t200);
	t171 = t149 * t156 * t192;
	t210 = (t150 * t167 + t152 * t171) / t142 ^ 2;
	t178 = 0.1e1 + t202;
	t126 = t178 * t156 * t140;
	t209 = t126 * t156;
	t208 = t128 * t144;
	t206 = t133 * t145;
	t203 = t149 * t151;
	t154 = t158 ^ 2;
	t201 = t149 * t154;
	t199 = t155 * t158;
	t193 = qJD(1) * t156;
	t120 = t122 * t201 + 0.1e1;
	t188 = 0.2e1 * (-t201 * t215 + (t154 * t155 * t190 - t171) * t122) / t120 ^ 2;
	t187 = -0.2e1 * t216;
	t186 = 0.2e1 * t215;
	t185 = t133 * t214;
	t184 = t122 * t199;
	t182 = t140 * t203;
	t177 = t155 * t188;
	t176 = t210 * t218;
	t175 = t210 * t220;
	t172 = t156 * t182;
	t170 = t178 * t158;
	t169 = t129 * t206 - t208;
	t166 = t155 * t191 + t173 * t158;
	t132 = -t145 * t195 + t204;
	t118 = 0.1e1 / t120;
	t115 = 0.1e1 / t117;
	t111 = (t217 * t155 * t136 - t137 * t172) * t158;
	t110 = -t136 * t195 + t205 + (t136 * t157 - t137 * t196) * t126;
	t109 = -t178 * t175 + (qJD(1) * t170 + t167 * t220) * t140;
	t105 = t187 + 0.2e1 * (-t112 * t115 * t129 + (-t115 * t214 - t129 * t216) * t133) * t133;
	t1 = [t151 * t158 * t176 + (qJD(2) * t170 - t193 * t200) * t140, t109, 0, 0, 0, 0; (t121 * t177 + (-t121 * t190 + (qJD(1) * t111 + t108) * t211) * t118) * t156 + (t122 * t177 * t111 + (-((t114 * t172 + t217 * t190 + t176) * t136 + (t175 * t203 - t212 + (t212 + (t218 - t219) * t191) * t140) * t137) * t184 + (-t122 * t190 + t155 * t186) * t111 + (-t121 + ((-t150 + t154) * t137 * t182 + t217 * t183) * t122) * t155 * qJD(1)) * t118) * t158, (t110 * t211 - t121 * t157) * t158 * t188 + ((-t121 * t193 + (-qJD(2) * t110 - t108) * t158 * t122) * t157 + (-t121 * t189 - (-t109 * t137 * t156 + t136 * t191 + t209 * t213 - t213 + (-qJD(2) * t136 - t137 * t192) * t126) * t184 + (t122 * t193 + t158 * t186) * t110 - ((t109 - t192) * t136 + ((0.1e1 - t209) * qJD(2) + (t126 - t156) * t114) * t137) * t122 * t194) * t155) * t118, 0, 0, 0, 0; 0.2e1 * (t128 * t168 + t132 * t207) * t216 + (0.2e1 * t132 * t185 - t174 * t128 * t197 + t166 * t208 + (-t174 * t133 * t198 + t132 * t112 + t113 * t168 - t166 * t206) * t129) * t115, t169 * t187 * t199 + (t169 * t157 * t189 + (-t169 * t193 + ((-t128 * t147 - 0.2e1 * t185) * t145 + (-t112 * t145 + (-t133 * t147 + t113) * t144) * t129) * t158) * t155) * t115, t105, 0, t105, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:22
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (2878->96), mult. (2949->204), div. (516->12), fcn. (3430->9), ass. (0->95)
	t161 = sin(qJ(2));
	t155 = t161 ^ 2;
	t163 = cos(qJ(2));
	t158 = 0.1e1 / t163 ^ 2;
	t208 = t155 * t158;
	t162 = sin(qJ(1));
	t226 = 0.2e1 * t162;
	t225 = t161 * t208;
	t152 = qJ(3) + pkin(11) + qJ(5) + qJ(6);
	t151 = cos(t152);
	t164 = cos(qJ(1));
	t200 = t163 * t164;
	t150 = sin(t152);
	t204 = t162 * t150;
	t140 = t151 * t200 + t204;
	t202 = t162 * t161;
	t145 = atan2(-t202, -t163);
	t144 = cos(t145);
	t143 = sin(t145);
	t189 = t143 * t202;
	t130 = -t144 * t163 - t189;
	t127 = 0.1e1 / t130;
	t134 = 0.1e1 / t140;
	t157 = 0.1e1 / t163;
	t128 = 0.1e1 / t130 ^ 2;
	t135 = 0.1e1 / t140 ^ 2;
	t224 = -0.2e1 * t161;
	t156 = t162 ^ 2;
	t149 = t156 * t208 + 0.1e1;
	t147 = 0.1e1 / t149;
	t223 = t147 - 0.1e1;
	t153 = qJD(3) + qJD(5) + qJD(6);
	t201 = t162 * t163;
	t174 = t150 * t201 + t151 * t164;
	t195 = qJD(2) * t164;
	t185 = t161 * t195;
	t118 = t174 * qJD(1) - t140 * t153 + t150 * t185;
	t203 = t162 * t151;
	t139 = t150 * t200 - t203;
	t133 = t139 ^ 2;
	t122 = t133 * t135 + 0.1e1;
	t213 = t135 * t139;
	t179 = -qJD(1) * t163 + t153;
	t180 = t153 * t163 - qJD(1);
	t210 = t150 * t164;
	t119 = -t180 * t210 + (t179 * t162 - t185) * t151;
	t220 = t119 * t134 * t135;
	t222 = (-t118 * t213 - t133 * t220) / t122 ^ 2;
	t198 = qJD(1) * t164;
	t186 = t161 * t198;
	t196 = qJD(2) * t163;
	t197 = qJD(2) * t162;
	t123 = (-(-t162 * t196 - t186) * t157 + t197 * t208) * t147;
	t211 = t144 * t161;
	t114 = (-t123 * t162 + qJD(2)) * t211 + (-t186 + (t123 - t197) * t163) * t143;
	t221 = t114 * t127 * t128;
	t219 = t123 * t143;
	t218 = t123 * t161;
	t217 = t128 * t161;
	t206 = t157 * t161;
	t173 = qJD(2) * (t157 * t225 + t206);
	t177 = t155 * t162 * t198;
	t216 = (t156 * t173 + t158 * t177) / t149 ^ 2;
	t184 = 0.1e1 + t208;
	t132 = t184 * t162 * t147;
	t215 = t132 * t162;
	t214 = t134 * t150;
	t212 = t139 * t151;
	t209 = t155 * t157;
	t160 = t164 ^ 2;
	t207 = t155 * t160;
	t205 = t161 * t164;
	t199 = qJD(1) * t162;
	t126 = t128 * t207 + 0.1e1;
	t194 = 0.2e1 * (-t207 * t221 + (t160 * t161 * t196 - t177) * t128) / t126 ^ 2;
	t193 = -0.2e1 * t222;
	t192 = 0.2e1 * t221;
	t191 = t139 * t220;
	t190 = t128 * t205;
	t188 = t147 * t209;
	t183 = t161 * t194;
	t182 = t216 * t224;
	t181 = t216 * t226;
	t178 = t162 * t188;
	t176 = t184 * t164;
	t175 = t135 * t212 - t214;
	t172 = t161 * t197 + t179 * t164;
	t138 = -t151 * t201 + t210;
	t124 = 0.1e1 / t126;
	t120 = 0.1e1 / t122;
	t117 = (t223 * t161 * t143 - t144 * t178) * t164;
	t116 = -t143 * t201 + t211 + (t143 * t163 - t144 * t202) * t132;
	t115 = -t184 * t181 + (qJD(1) * t176 + t173 * t226) * t147;
	t111 = t193 + 0.2e1 * (-t118 * t120 * t135 + (-t120 * t220 - t135 * t222) * t139) * t139;
	t1 = [t157 * t164 * t182 + (qJD(2) * t176 - t199 * t206) * t147, t115, 0, 0, 0, 0; (t127 * t183 + (-t127 * t196 + (qJD(1) * t117 + t114) * t217) * t124) * t162 + (t128 * t183 * t117 + (-((t123 * t178 + t223 * t196 + t182) * t143 + (t181 * t209 - t218 + (t218 + (t224 - t225) * t197) * t147) * t144) * t190 + (-t128 * t196 + t161 * t192) * t117 + (-t127 + ((-t156 + t160) * t144 * t188 + t223 * t189) * t128) * t161 * qJD(1)) * t124) * t164, (t116 * t217 - t127 * t163) * t164 * t194 + ((-t127 * t199 + (-qJD(2) * t116 - t114) * t164 * t128) * t163 + (-t127 * t195 - (-t115 * t144 * t162 + t143 * t197 + t215 * t219 - t219 + (-qJD(2) * t143 - t144 * t198) * t132) * t190 + (t128 * t199 + t164 * t192) * t116 - ((t115 - t198) * t143 + ((0.1e1 - t215) * qJD(2) + (t132 - t162) * t123) * t144) * t128 * t200) * t161) * t124, 0, 0, 0, 0; 0.2e1 * (t134 * t174 + t138 * t213) * t222 + (0.2e1 * t138 * t191 - t180 * t134 * t203 + t172 * t214 + (-t180 * t139 * t204 + t138 * t118 + t119 * t174 - t172 * t212) * t135) * t120, t175 * t193 * t205 + (t175 * t163 * t195 + (-t175 * t199 + ((-t134 * t153 - 0.2e1 * t191) * t151 + (-t118 * t151 + (-t139 * t153 + t119) * t150) * t135) * t164) * t161) * t120, t111, 0, t111, t111;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end