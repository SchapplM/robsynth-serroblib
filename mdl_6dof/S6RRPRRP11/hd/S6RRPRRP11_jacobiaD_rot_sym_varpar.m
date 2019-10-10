% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP11
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
%   Wie in S6RRPRRP11_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP11_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP11_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:44
	% DurationCPUTime: 0.78s
	% Computational Cost: add. (774->70), mult. (1839->155), div. (436->14), fcn. (2165->7), ass. (0->75)
	t87 = cos(qJ(1));
	t116 = qJD(1) * t87;
	t85 = sin(qJ(1));
	t138 = 0.2e1 * t85;
	t74 = t85 ^ 2;
	t75 = 0.1e1 / t85;
	t83 = t87 ^ 2;
	t136 = (0.1e1 + 0.1e1 / t74 * t83) * t75 * t116;
	t84 = sin(qJ(2));
	t118 = t85 * t84;
	t86 = cos(qJ(2));
	t65 = atan2(-t118, -t86);
	t63 = sin(t65);
	t108 = t63 * t118;
	t64 = cos(t65);
	t59 = -t64 * t86 - t108;
	t56 = 0.1e1 / t59;
	t79 = 0.1e1 / t86;
	t57 = 0.1e1 / t59 ^ 2;
	t135 = -0.2e1 * t84;
	t73 = t84 ^ 2;
	t80 = 0.1e1 / t86 ^ 2;
	t123 = t73 * t80;
	t70 = t74 * t123 + 0.1e1;
	t66 = 0.1e1 / t70;
	t134 = t66 - 0.1e1;
	t106 = t84 * t116;
	t115 = qJD(2) * t85;
	t125 = t64 * t84;
	t114 = qJD(2) * t86;
	t52 = (-(-t85 * t114 - t106) * t79 + t115 * t123) * t66;
	t48 = (-t52 * t85 + qJD(2)) * t125 + (-t106 + (t52 - t115) * t86) * t63;
	t133 = t48 * t56 * t57;
	t132 = t52 * t63;
	t131 = t52 * t84;
	t130 = t57 * t84;
	t129 = t57 * t87;
	t120 = t79 * t84;
	t72 = t84 * t73;
	t78 = t86 ^ 2;
	t94 = qJD(2) * (t72 * t79 / t78 + t120);
	t99 = t73 * t85 * t116;
	t128 = (t74 * t94 + t80 * t99) / t70 ^ 2;
	t105 = 0.1e1 + t123;
	t62 = t105 * t85 * t66;
	t127 = t62 * t85;
	t126 = t63 * t86;
	t124 = t73 * t79;
	t122 = t73 * t83;
	t76 = 0.1e1 / t85 ^ 2;
	t121 = t76 * t83;
	t119 = t84 * t87;
	t117 = qJD(1) * t85;
	t113 = qJD(2) * t87;
	t55 = t57 * t122 + 0.1e1;
	t98 = t83 * t84 * t114;
	t112 = 0.2e1 * (-t122 * t133 + (t98 - t99) * t57) / t55 ^ 2;
	t111 = 0.2e1 * t133;
	t71 = t78 * t121 + 0.1e1;
	t110 = 0.2e1 * (-t78 * t136 - t76 * t98) / t71 ^ 2;
	t109 = t57 * t119;
	t107 = t66 * t124;
	t104 = 0.1e1 + t121;
	t103 = t84 * t112;
	t102 = t128 * t135;
	t101 = t128 * t138;
	t100 = t85 * t107;
	t97 = t105 * t87;
	t96 = t104 * t84;
	t68 = 0.1e1 / t71;
	t53 = 0.1e1 / t55;
	t51 = (t134 * t84 * t63 - t64 * t100) * t87;
	t50 = -t85 * t126 + t125 + (-t64 * t118 + t126) * t62;
	t49 = -t105 * t101 + (qJD(1) * t97 + t94 * t138) * t66;
	t1 = [t87 * t79 * t102 + (qJD(2) * t97 - t117 * t120) * t66, t49, 0, 0, 0, 0; (t56 * t103 + (-t56 * t114 + (qJD(1) * t51 + t48) * t130) * t53) * t85 + (t57 * t103 * t51 + (-((t52 * t100 + t134 * t114 + t102) * t63 + (t101 * t124 - t131 + (t131 + (-t72 * t80 + t135) * t115) * t66) * t64) * t109 + (t84 * t111 - t57 * t114) * t51 + (-t56 + ((-t74 + t83) * t64 * t107 + t134 * t108) * t57) * t84 * qJD(1)) * t53) * t87, (t50 * t130 - t56 * t86) * t87 * t112 + ((-t56 * t117 + (-qJD(2) * t50 - t48) * t129) * t86 + (-t56 * t113 - (-t49 * t64 * t85 + t63 * t115 + t127 * t132 - t132 + (-qJD(2) * t63 - t116 * t64) * t62) * t109 + (t87 * t111 + t57 * t117) * t50 - ((t49 - t116) * t63 + ((0.1e1 - t127) * qJD(2) + (t62 - t85) * t52) * t64) * t86 * t129) * t84) * t53, 0, 0, 0, 0; t104 * t86 * t110 + (qJD(2) * t96 + 0.2e1 * t86 * t136) * t68, t75 * t110 * t119 + (-t75 * t86 * t113 + qJD(1) * t96) * t68, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:44
	% DurationCPUTime: 0.96s
	% Computational Cost: add. (813->90), mult. (2519->205), div. (480->12), fcn. (2968->9), ass. (0->93)
	t131 = cos(qJ(2));
	t129 = sin(qJ(1));
	t123 = t129 ^ 2;
	t128 = sin(qJ(2));
	t121 = 0.1e1 / t128 ^ 2;
	t125 = t131 ^ 2;
	t179 = t121 * t125;
	t118 = t123 * t179 + 0.1e1;
	t116 = 0.1e1 / t118;
	t120 = 0.1e1 / t128;
	t169 = qJD(2) * t129;
	t157 = t121 * t169;
	t132 = cos(qJ(1));
	t171 = qJD(1) * t132;
	t158 = t131 * t171;
	t90 = ((t128 * t169 - t158) * t120 + t125 * t157) * t116;
	t193 = t131 * t90;
	t148 = -t90 + t169;
	t146 = qJD(1) * t128 + qJD(4);
	t167 = qJD(2) * t132;
	t192 = t146 * t129 - t131 * t167;
	t176 = t129 * t131;
	t115 = atan2(-t176, t128);
	t114 = cos(t115);
	t113 = sin(t115);
	t161 = t113 * t176;
	t99 = t114 * t128 - t161;
	t96 = 0.1e1 / t99;
	t127 = sin(qJ(4));
	t175 = t132 * t127;
	t130 = cos(qJ(4));
	t177 = t129 * t130;
	t110 = t128 * t175 + t177;
	t106 = 0.1e1 / t110;
	t107 = 0.1e1 / t110 ^ 2;
	t97 = 0.1e1 / t99 ^ 2;
	t191 = t116 - 0.1e1;
	t149 = -t129 * t90 + qJD(2);
	t181 = t114 * t131;
	t85 = t149 * t181 + (t148 * t128 - t158) * t113;
	t190 = t85 * t96 * t97;
	t126 = t132 ^ 2;
	t95 = t126 * t125 * t97 + 0.1e1;
	t93 = 0.1e1 / t95;
	t189 = t93 * t97;
	t188 = t96 * t93;
	t174 = t132 * t130;
	t178 = t129 * t127;
	t109 = -t128 * t174 + t178;
	t105 = t109 ^ 2;
	t104 = t105 * t107 + 0.1e1;
	t184 = t107 * t109;
	t147 = qJD(4) * t128 + qJD(1);
	t143 = t147 * t132;
	t92 = -t192 * t127 + t130 * t143;
	t186 = t106 * t107 * t92;
	t91 = t127 * t143 + t192 * t130;
	t187 = 0.1e1 / t104 ^ 2 * (-t105 * t186 + t91 * t184);
	t185 = t132 * t97;
	t183 = t109 * t127;
	t182 = t113 * t129;
	t180 = t120 * t125;
	t173 = qJD(1) * t129;
	t172 = qJD(1) * t131;
	t170 = qJD(2) * t128;
	t168 = qJD(2) * t131;
	t159 = t129 * t171;
	t162 = t97 * t170;
	t166 = 0.2e1 * (-t126 * t131 * t162 + (-t126 * t190 - t97 * t159) * t125) / t95 ^ 2;
	t165 = 0.2e1 * t190;
	t164 = 0.2e1 * t187;
	t124 = t131 * t125;
	t141 = qJD(2) * (-t121 * t124 - t131) * t120;
	t163 = 0.2e1 * (t123 * t141 + t159 * t179) / t118 ^ 2;
	t160 = t116 * t180;
	t155 = t96 * t166;
	t154 = t97 * t166;
	t153 = 0.1e1 + t179;
	t152 = 0.2e1 * t109 * t186;
	t151 = t129 * t163;
	t150 = t131 * t163;
	t145 = t129 * t160;
	t144 = t153 * t132;
	t142 = t106 * t130 + t107 * t183;
	t140 = t142 * t132;
	t112 = -t128 * t178 + t174;
	t111 = t128 * t177 + t175;
	t103 = t153 * t129 * t116;
	t101 = 0.1e1 / t104;
	t89 = (t191 * t131 * t113 + t114 * t145) * t132;
	t88 = t128 * t182 + t181 + (-t113 * t128 - t114 * t176) * t103;
	t86 = -t153 * t151 + (qJD(1) * t144 + 0.2e1 * t129 * t141) * t116;
	t1 = [t132 * t120 * t150 + (t120 * t129 * t172 + qJD(2) * t144) * t116, t86, 0, 0, 0, 0; (t170 * t188 + (t155 + (qJD(1) * t89 + t85) * t189) * t131) * t129 + (t89 * t154 * t131 + (t89 * t162 + (t89 * t165 + ((t90 * t145 + t191 * t170 + t150) * t113 + (t151 * t180 + t193 + (t124 * t157 - (t90 - 0.2e1 * t169) * t131) * t116) * t114) * t185) * t131 + (-t96 + (-(-t123 + t126) * t114 * t160 + t191 * t161) * t97) * t172) * t93) * t132, (t173 * t188 + (t155 + (qJD(2) * t88 + t85) * t189) * t132) * t128 + (t88 * t132 * t154 + (-t96 * t167 + (t132 * t165 + t97 * t173) * t88 + (-t103 * t182 * t193 + (-(-t103 * t171 - t129 * t86) * t131 - (t148 * t103 - t149) * t128) * t114 + (-(-qJD(2) * t103 + t148) * t131 - (-t86 + t171) * t128) * t113) * t185) * t93) * t131, 0, 0, 0, 0; (-t106 * t111 + t112 * t184) * t164 + (t112 * t152 + (-t111 * t92 - t112 * t91 + t147 * t109 * t177 - (-t129 * t168 - t146 * t132) * t183) * t107 + (t146 * t174 + (-t147 * t127 + t130 * t168) * t129) * t106) * t101, t131 * t140 * t164 + (t140 * t170 + (t142 * t173 + ((qJD(4) * t106 + t152) * t127 + (-t127 * t91 + (-qJD(4) * t109 + t92) * t130) * t107) * t132) * t131) * t101, 0, -0.2e1 * t187 + 0.2e1 * (t91 * t107 * t101 + (-t101 * t186 - t107 * t187) * t109) * t109, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:44
	% DurationCPUTime: 1.04s
	% Computational Cost: add. (1381->94), mult. (2734->208), div. (498->12), fcn. (3199->9), ass. (0->96)
	t152 = sin(qJ(2));
	t145 = 0.1e1 / t152 ^ 2;
	t154 = cos(qJ(2));
	t149 = t154 ^ 2;
	t200 = t145 * t149;
	t218 = t154 * t200;
	t153 = sin(qJ(1));
	t175 = 0.1e1 + t200;
	t217 = t153 * t175;
	t143 = qJD(4) + qJD(5);
	t171 = qJD(1) * t152 + t143;
	t155 = cos(qJ(1));
	t187 = qJD(2) * t155;
	t216 = t153 * t171 - t154 * t187;
	t188 = qJD(2) * t154;
	t215 = t153 * t188 + t155 * t171;
	t194 = t153 * t154;
	t136 = atan2(-t194, t152);
	t135 = cos(t136);
	t134 = sin(t136);
	t180 = t134 * t194;
	t123 = t135 * t152 - t180;
	t120 = 0.1e1 / t123;
	t151 = qJ(4) + qJ(5);
	t141 = sin(t151);
	t142 = cos(t151);
	t195 = t153 * t142;
	t196 = t152 * t155;
	t131 = t141 * t196 + t195;
	t127 = 0.1e1 / t131;
	t144 = 0.1e1 / t152;
	t121 = 0.1e1 / t123 ^ 2;
	t128 = 0.1e1 / t131 ^ 2;
	t147 = t153 ^ 2;
	t140 = t147 * t200 + 0.1e1;
	t137 = 0.1e1 / t140;
	t214 = t137 - 0.1e1;
	t172 = t143 * t152 + qJD(1);
	t166 = t172 * t155;
	t111 = t141 * t166 + t142 * t216;
	t130 = t141 * t153 - t142 * t196;
	t126 = t130 ^ 2;
	t119 = t126 * t128 + 0.1e1;
	t203 = t128 * t130;
	t112 = -t141 * t216 + t142 * t166;
	t211 = t112 * t127 * t128;
	t213 = (t111 * t203 - t126 * t211) / t119 ^ 2;
	t191 = qJD(1) * t155;
	t178 = t154 * t191;
	t189 = qJD(2) * t153;
	t113 = ((t152 * t189 - t178) * t144 + t189 * t200) * t137;
	t201 = t135 * t154;
	t107 = (-t113 * t153 + qJD(2)) * t201 + (-t178 + (-t113 + t189) * t152) * t134;
	t212 = t107 * t120 * t121;
	t210 = t113 * t134;
	t209 = t113 * t154;
	t208 = t121 * t154;
	t207 = t121 * t155;
	t164 = qJD(2) * (-t154 - t218) * t144;
	t198 = t149 * t153;
	t169 = t191 * t198;
	t206 = (t145 * t169 + t147 * t164) / t140 ^ 2;
	t125 = t137 * t217;
	t205 = t125 * t153;
	t204 = t127 * t142;
	t202 = t130 * t141;
	t150 = t155 ^ 2;
	t199 = t149 * t150;
	t197 = t152 * t153;
	t193 = qJD(1) * t153;
	t192 = qJD(1) * t154;
	t190 = qJD(2) * t152;
	t116 = t121 * t199 + 0.1e1;
	t186 = 0.2e1 * (-t199 * t212 + (-t150 * t152 * t188 - t169) * t121) / t116 ^ 2;
	t185 = 0.2e1 * t213;
	t184 = 0.2e1 * t212;
	t183 = -0.2e1 * t206;
	t182 = t154 * t207;
	t181 = t154 * t206;
	t179 = t144 * t198;
	t174 = t154 * t186;
	t173 = 0.2e1 * t130 * t211;
	t170 = t137 * t179;
	t168 = t175 * t155;
	t167 = t172 * t153;
	t165 = t128 * t202 + t204;
	t163 = t165 * t155;
	t133 = -t141 * t197 + t142 * t155;
	t132 = t141 * t155 + t152 * t195;
	t117 = 0.1e1 / t119;
	t114 = 0.1e1 / t116;
	t110 = (t134 * t154 * t214 + t135 * t170) * t155;
	t109 = t134 * t197 + t201 + (-t134 * t152 - t135 * t194) * t125;
	t108 = t183 * t217 + (qJD(1) * t168 + 0.2e1 * t153 * t164) * t137;
	t104 = -0.2e1 * t213 + 0.2e1 * (t111 * t117 * t128 + (-t117 * t211 - t128 * t213) * t130) * t130;
	t1 = [0.2e1 * t144 * t155 * t181 + (t144 * t153 * t192 + qJD(2) * t168) * t137, t108, 0, 0, 0, 0; (t120 * t174 + (t120 * t190 + (qJD(1) * t110 + t107) * t208) * t114) * t153 + (t121 * t174 * t110 + (-((-t113 * t170 - t190 * t214 - 0.2e1 * t181) * t134 + (t179 * t183 - t209 + (t209 + (-0.2e1 * t154 - t218) * t189) * t137) * t135) * t182 + (t121 * t190 + t154 * t184) * t110 + (-t120 + ((t147 - t150) * t149 * t144 * t137 * t135 + t214 * t180) * t121) * t192) * t114) * t155, (t109 * t208 + t120 * t152) * t155 * t186 + ((t120 * t193 + (qJD(2) * t109 + t107) * t207) * t152 + (-t120 * t187 - (-t108 * t135 * t153 + t134 * t189 + t205 * t210 - t210 + (-qJD(2) * t134 - t135 * t191) * t125) * t182 + (t121 * t193 + t155 * t184) * t109 - ((-t108 + t191) * t134 + ((-0.1e1 + t205) * qJD(2) + (-t125 + t153) * t113) * t135) * t121 * t196) * t154) * t114, 0, 0, 0, 0; (-t127 * t132 + t133 * t203) * t185 + (t133 * t173 - t127 * t141 * t167 + t215 * t204 + (t130 * t142 * t167 - t133 * t111 - t132 * t112 + t202 * t215) * t128) * t117, t154 * t163 * t185 + (t163 * t190 + (t165 * t193 + ((t127 * t143 + t173) * t141 + (-t111 * t141 + (-t130 * t143 + t112) * t142) * t128) * t155) * t154) * t117, 0, t104, t104, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:44
	% DurationCPUTime: 1.04s
	% Computational Cost: add. (1381->94), mult. (2734->208), div. (498->12), fcn. (3199->9), ass. (0->96)
	t152 = sin(qJ(2));
	t145 = 0.1e1 / t152 ^ 2;
	t154 = cos(qJ(2));
	t149 = t154 ^ 2;
	t200 = t145 * t149;
	t218 = t154 * t200;
	t153 = sin(qJ(1));
	t175 = 0.1e1 + t200;
	t217 = t153 * t175;
	t143 = qJD(4) + qJD(5);
	t171 = qJD(1) * t152 + t143;
	t155 = cos(qJ(1));
	t187 = qJD(2) * t155;
	t216 = t153 * t171 - t154 * t187;
	t188 = qJD(2) * t154;
	t215 = t153 * t188 + t155 * t171;
	t194 = t153 * t154;
	t136 = atan2(-t194, t152);
	t135 = cos(t136);
	t134 = sin(t136);
	t180 = t134 * t194;
	t123 = t135 * t152 - t180;
	t120 = 0.1e1 / t123;
	t151 = qJ(4) + qJ(5);
	t141 = sin(t151);
	t142 = cos(t151);
	t195 = t153 * t142;
	t196 = t152 * t155;
	t131 = t141 * t196 + t195;
	t127 = 0.1e1 / t131;
	t144 = 0.1e1 / t152;
	t121 = 0.1e1 / t123 ^ 2;
	t128 = 0.1e1 / t131 ^ 2;
	t147 = t153 ^ 2;
	t140 = t147 * t200 + 0.1e1;
	t137 = 0.1e1 / t140;
	t214 = t137 - 0.1e1;
	t172 = t143 * t152 + qJD(1);
	t166 = t172 * t155;
	t111 = t141 * t166 + t142 * t216;
	t130 = t141 * t153 - t142 * t196;
	t126 = t130 ^ 2;
	t119 = t126 * t128 + 0.1e1;
	t203 = t128 * t130;
	t112 = -t141 * t216 + t142 * t166;
	t211 = t112 * t127 * t128;
	t213 = (t111 * t203 - t126 * t211) / t119 ^ 2;
	t191 = qJD(1) * t155;
	t178 = t154 * t191;
	t189 = qJD(2) * t153;
	t113 = ((t152 * t189 - t178) * t144 + t189 * t200) * t137;
	t201 = t135 * t154;
	t107 = (-t113 * t153 + qJD(2)) * t201 + (-t178 + (-t113 + t189) * t152) * t134;
	t212 = t107 * t120 * t121;
	t210 = t113 * t134;
	t209 = t113 * t154;
	t208 = t121 * t154;
	t207 = t121 * t155;
	t164 = qJD(2) * (-t154 - t218) * t144;
	t198 = t149 * t153;
	t169 = t191 * t198;
	t206 = (t145 * t169 + t147 * t164) / t140 ^ 2;
	t125 = t137 * t217;
	t205 = t125 * t153;
	t204 = t127 * t142;
	t202 = t130 * t141;
	t150 = t155 ^ 2;
	t199 = t149 * t150;
	t197 = t152 * t153;
	t193 = qJD(1) * t153;
	t192 = qJD(1) * t154;
	t190 = qJD(2) * t152;
	t116 = t121 * t199 + 0.1e1;
	t186 = 0.2e1 * (-t199 * t212 + (-t150 * t152 * t188 - t169) * t121) / t116 ^ 2;
	t185 = 0.2e1 * t213;
	t184 = 0.2e1 * t212;
	t183 = -0.2e1 * t206;
	t182 = t154 * t207;
	t181 = t154 * t206;
	t179 = t144 * t198;
	t174 = t154 * t186;
	t173 = 0.2e1 * t130 * t211;
	t170 = t137 * t179;
	t168 = t175 * t155;
	t167 = t172 * t153;
	t165 = t128 * t202 + t204;
	t163 = t165 * t155;
	t133 = -t141 * t197 + t142 * t155;
	t132 = t141 * t155 + t152 * t195;
	t117 = 0.1e1 / t119;
	t114 = 0.1e1 / t116;
	t110 = (t134 * t154 * t214 + t135 * t170) * t155;
	t109 = t134 * t197 + t201 + (-t134 * t152 - t135 * t194) * t125;
	t108 = t183 * t217 + (qJD(1) * t168 + 0.2e1 * t153 * t164) * t137;
	t104 = -0.2e1 * t213 + 0.2e1 * (t111 * t117 * t128 + (-t117 * t211 - t128 * t213) * t130) * t130;
	t1 = [0.2e1 * t144 * t155 * t181 + (t144 * t153 * t192 + qJD(2) * t168) * t137, t108, 0, 0, 0, 0; (t120 * t174 + (t120 * t190 + (qJD(1) * t110 + t107) * t208) * t114) * t153 + (t121 * t174 * t110 + (-((-t113 * t170 - t190 * t214 - 0.2e1 * t181) * t134 + (t179 * t183 - t209 + (t209 + (-0.2e1 * t154 - t218) * t189) * t137) * t135) * t182 + (t121 * t190 + t154 * t184) * t110 + (-t120 + ((t147 - t150) * t149 * t144 * t137 * t135 + t214 * t180) * t121) * t192) * t114) * t155, (t109 * t208 + t120 * t152) * t155 * t186 + ((t120 * t193 + (qJD(2) * t109 + t107) * t207) * t152 + (-t120 * t187 - (-t108 * t135 * t153 + t134 * t189 + t205 * t210 - t210 + (-qJD(2) * t134 - t135 * t191) * t125) * t182 + (t121 * t193 + t155 * t184) * t109 - ((-t108 + t191) * t134 + ((-0.1e1 + t205) * qJD(2) + (-t125 + t153) * t113) * t135) * t121 * t196) * t154) * t114, 0, 0, 0, 0; (-t127 * t132 + t133 * t203) * t185 + (t133 * t173 - t127 * t141 * t167 + t215 * t204 + (t130 * t142 * t167 - t133 * t111 - t132 * t112 + t202 * t215) * t128) * t117, t154 * t163 * t185 + (t163 * t190 + (t165 * t193 + ((t127 * t143 + t173) * t141 + (-t111 * t141 + (-t130 * t143 + t112) * t142) * t128) * t155) * t154) * t117, 0, t104, t104, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end