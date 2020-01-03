% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPR9
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
%   Wie in S5RRRPR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRPR9_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR9_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:25:42
	% EndTime: 2019-12-31 21:25:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:25:42
	% EndTime: 2019-12-31 21:25:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:25:42
	% EndTime: 2019-12-31 21:25:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:25:42
	% EndTime: 2019-12-31 21:25:43
	% DurationCPUTime: 0.72s
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
	t139 = t102 * t175 - t177;
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
	t188 = (t123 * t125 * t155 + (-t123 * t187 - t92 * t151) * t118) / t90 ^ 2;
	t176 = t102 * t106;
	t143 = -qJD(1) * t128 + qJD(3);
	t144 = qJD(3) * t128 - qJD(1);
	t160 = qJD(2) * t129;
	t148 = t125 * t160;
	t87 = -t144 * t166 + (t143 * t126 - t148) * t127;
	t183 = t101 * t102 * t87;
	t167 = t126 * t128;
	t138 = t124 * t167 + t165;
	t86 = t138 * qJD(1) - t107 * qJD(3) + t124 * t148;
	t186 = (-t100 * t183 - t86 * t176) / t99 ^ 2;
	t88 = 0.1e1 / t90;
	t185 = t88 * t91;
	t184 = t88 * t92;
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
	t84 = (t189 * t125 * t108 - t109 * t142) * t129;
	t83 = -t192 * t174 + (-t126 + t98) * t128 * t108;
	t81 = t190 * t193 + (qJD(1) * t141 + 0.2e1 * t126 * t137) * t111;
	t1 = [t145 * t171 + (qJD(2) * t141 - t164 * t172) * t111, t81, 0, 0, 0; (-t161 * t185 + (0.2e1 * t157 + (qJD(1) * t84 + t80) * t184) * t125) * t126 + (t84 * t147 * t125 + (-t84 * t155 + (t84 * t159 + ((0.2e1 * t125 * t181 - t85 * t142 - t189 * t161) * t108 + (t118 * t126 * t145 + t125 * t85 + (t117 * t149 - (t85 - 0.2e1 * t162) * t125) * t111) * t109) * t179) * t125 + (-t91 + (-(t119 - t123) * t109 * t153 + t189 * t154) * t92) * t125 * qJD(1)) * t88) * t129, (-t164 * t185 + (-0.2e1 * t157 + (-qJD(2) * t83 - t80) * t184) * t129) * t128 + (t83 * t129 * t147 + (-t91 * t160 - ((-t126 * t81 - t163 * t98) * t109 + (t192 * t85 + t162 - t178) * t108) * t92 * t171 + (t129 * t159 + t92 * t164) * t83 - ((t81 - t163) * t108 + (t85 * t98 + qJD(2) + (-t85 - t178) * t126) * t109) * t128 * t179) * t88) * t125, 0, 0, 0; 0.2e1 * (t101 * t138 + t105 * t176) * t186 + (0.2e1 * t105 * t156 - t144 * t101 * t168 + (t125 * t162 + t140) * t177 + (t138 * t87 + t105 * t86 - t140 * t175 - (qJD(2) * t125 * t127 + t144 * t124) * t106 * t126) * t102) * t96, t128 * t160 * t191 + (-t164 * t191 + (t139 * t158 + ((-qJD(3) * t101 - 0.2e1 * t156) * t127 + (-t127 * t86 + (-qJD(3) * t106 + t87) * t124) * t102) * t96) * t129) * t125, t158 + 0.2e1 * (-t86 * t102 * t96 + (-t102 * t186 - t96 * t183) * t106) * t106, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:25:42
	% EndTime: 2019-12-31 21:25:43
	% DurationCPUTime: 0.72s
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
	t126 = qJ(3) + pkin(9);
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
	t1 = [t156 * t180 + (qJD(2) * t149 - t174 * t181) * t119, t89, 0, 0, 0; (-t101 * t164 + (0.2e1 * t166 + (qJD(1) * t91 + t88) * t194) * t134) * t135 + ((-t91 * t164 + (t91 * t169 + ((0.2e1 * t134 * t190 - t94 * t151 - t197 * t171) * t115 + (t156 * t182 + t134 * t94 + (t127 * t159 - (t94 - 0.2e1 * t172) * t134) * t119) * t116) * t95 * t137) * t134) * t102 + (t91 * t167 + (-t101 + ((-t129 + t133) * t116 * t162 + t197 * t163) * t102) * qJD(1)) * t134 * t95) * t137, (-t101 * t95 * t174 + (-0.2e1 * t166 + (-qJD(2) * t90 - t88) * t194) * t137) * t136 + (t90 * t137 * t102 * t169 + ((-qJD(2) * t101 + t90 * t167) * t137 + (t90 * t174 + (-(-t106 * t173 - t135 * t89) * t116 - ((t106 * t135 - 0.1e1) * t94 + (-t106 + t135) * qJD(2)) * t115) * t180) * t102) * t95 - ((t89 - t173) * t115 + (t154 * t106 + t155) * t116) * t175 * t194) * t134, 0, 0, 0; 0.2e1 * (t108 * t147 + t112 * t188) * t195 + (0.2e1 * t112 * t165 - t153 * t108 * t178 + t145 * t189 + (-t153 * t113 * t179 + t112 * t92 - t145 * t187 + t147 * t93) * t109) * t98, t136 * t170 * t199 + (-t174 * t199 + (t148 * t168 + ((-qJD(3) * t108 - 0.2e1 * t165) * t125 + (-t125 * t92 + (-qJD(3) * t113 + t93) * t124) * t109) * t98) * t137) * t134, t168 + 0.2e1 * (-t109 * t92 * t98 + (-t109 * t195 - t98 * t192) * t113) * t113, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:25:42
	% EndTime: 2019-12-31 21:25:43
	% DurationCPUTime: 0.76s
	% Computational Cost: add. (2009->96), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->95)
	t154 = sin(qJ(2));
	t148 = t154 ^ 2;
	t156 = cos(qJ(2));
	t151 = 0.1e1 / t156 ^ 2;
	t201 = t148 * t151;
	t155 = sin(qJ(1));
	t219 = 0.2e1 * t155;
	t218 = t154 * t201;
	t145 = qJ(3) + pkin(9) + qJ(5);
	t144 = cos(t145);
	t157 = cos(qJ(1));
	t193 = t156 * t157;
	t143 = sin(t145);
	t197 = t155 * t143;
	t133 = t144 * t193 + t197;
	t195 = t155 * t154;
	t138 = atan2(-t195, -t156);
	t136 = cos(t138);
	t135 = sin(t138);
	t182 = t135 * t195;
	t123 = -t136 * t156 - t182;
	t120 = 0.1e1 / t123;
	t127 = 0.1e1 / t133;
	t150 = 0.1e1 / t156;
	t121 = 0.1e1 / t123 ^ 2;
	t128 = 0.1e1 / t133 ^ 2;
	t217 = -0.2e1 * t154;
	t149 = t155 ^ 2;
	t141 = t149 * t201 + 0.1e1;
	t139 = 0.1e1 / t141;
	t216 = t139 - 0.1e1;
	t146 = qJD(3) + qJD(5);
	t194 = t155 * t156;
	t167 = t143 * t194 + t144 * t157;
	t188 = qJD(2) * t157;
	t178 = t154 * t188;
	t111 = t167 * qJD(1) - t133 * t146 + t143 * t178;
	t196 = t155 * t144;
	t132 = t143 * t193 - t196;
	t126 = t132 ^ 2;
	t116 = t126 * t128 + 0.1e1;
	t206 = t128 * t132;
	t172 = -qJD(1) * t156 + t146;
	t173 = t146 * t156 - qJD(1);
	t203 = t143 * t157;
	t112 = -t173 * t203 + (t172 * t155 - t178) * t144;
	t213 = t112 * t127 * t128;
	t215 = (-t111 * t206 - t126 * t213) / t116 ^ 2;
	t191 = qJD(1) * t157;
	t179 = t154 * t191;
	t189 = qJD(2) * t156;
	t190 = qJD(2) * t155;
	t113 = (-(-t155 * t189 - t179) * t150 + t190 * t201) * t139;
	t204 = t136 * t154;
	t107 = (-t113 * t155 + qJD(2)) * t204 + (-t179 + (t113 - t190) * t156) * t135;
	t214 = t107 * t120 * t121;
	t212 = t113 * t135;
	t211 = t113 * t154;
	t210 = t121 * t154;
	t199 = t150 * t154;
	t166 = qJD(2) * (t150 * t218 + t199);
	t170 = t148 * t155 * t191;
	t209 = (t149 * t166 + t151 * t170) / t141 ^ 2;
	t177 = 0.1e1 + t201;
	t125 = t177 * t155 * t139;
	t208 = t125 * t155;
	t207 = t127 * t143;
	t205 = t132 * t144;
	t202 = t148 * t150;
	t153 = t157 ^ 2;
	t200 = t148 * t153;
	t198 = t154 * t157;
	t192 = qJD(1) * t155;
	t119 = t121 * t200 + 0.1e1;
	t187 = 0.2e1 * (-t200 * t214 + (t153 * t154 * t189 - t170) * t121) / t119 ^ 2;
	t186 = -0.2e1 * t215;
	t185 = 0.2e1 * t214;
	t184 = t132 * t213;
	t183 = t121 * t198;
	t181 = t139 * t202;
	t176 = t154 * t187;
	t175 = t209 * t217;
	t174 = t209 * t219;
	t171 = t155 * t181;
	t169 = t177 * t157;
	t168 = t128 * t205 - t207;
	t165 = t154 * t190 + t172 * t157;
	t131 = -t144 * t194 + t203;
	t117 = 0.1e1 / t119;
	t114 = 0.1e1 / t116;
	t110 = (t216 * t154 * t135 - t136 * t171) * t157;
	t109 = -t135 * t194 + t204 + (t135 * t156 - t136 * t195) * t125;
	t108 = -t177 * t174 + (qJD(1) * t169 + t166 * t219) * t139;
	t104 = t186 + 0.2e1 * (-t111 * t114 * t128 + (-t114 * t213 - t128 * t215) * t132) * t132;
	t1 = [t150 * t157 * t175 + (qJD(2) * t169 - t192 * t199) * t139, t108, 0, 0, 0; (t120 * t176 + (-t120 * t189 + (qJD(1) * t110 + t107) * t210) * t117) * t155 + (t121 * t176 * t110 + (-((t113 * t171 + t216 * t189 + t175) * t135 + (t174 * t202 - t211 + (t211 + (t217 - t218) * t190) * t139) * t136) * t183 + (-t121 * t189 + t154 * t185) * t110 + (-t120 + ((-t149 + t153) * t136 * t181 + t216 * t182) * t121) * t154 * qJD(1)) * t117) * t157, (t109 * t210 - t120 * t156) * t157 * t187 + ((-t120 * t192 + (-qJD(2) * t109 - t107) * t157 * t121) * t156 + (-t120 * t188 - (-t108 * t136 * t155 + t135 * t190 + t208 * t212 - t212 + (-qJD(2) * t135 - t136 * t191) * t125) * t183 + (t121 * t192 + t157 * t185) * t109 - ((t108 - t191) * t135 + ((0.1e1 - t208) * qJD(2) + (t125 - t155) * t113) * t136) * t121 * t193) * t154) * t117, 0, 0, 0; 0.2e1 * (t127 * t167 + t131 * t206) * t215 + (0.2e1 * t131 * t184 - t173 * t127 * t196 + t165 * t207 + (-t173 * t132 * t197 + t131 * t111 + t112 * t167 - t165 * t205) * t128) * t114, t168 * t186 * t198 + (t168 * t156 * t188 + (-t168 * t192 + ((-t127 * t146 - 0.2e1 * t184) * t144 + (-t111 * t144 + (-t132 * t146 + t112) * t143) * t128) * t157) * t154) * t114, t104, 0, t104;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end