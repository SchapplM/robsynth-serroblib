% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR2
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
%   Wie in S6RRPPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:19
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:54
	% EndTime: 2019-10-10 09:19:54
	% DurationCPUTime: 0.77s
	% Computational Cost: add. (1883->70), mult. (1839->155), div. (436->14), fcn. (2165->7), ass. (0->75)
	t95 = sin(qJ(1));
	t145 = 0.2e1 * t95;
	t90 = 0.1e1 / t95;
	t96 = cos(qJ(1));
	t129 = t90 * t96;
	t123 = qJD(2) * t95;
	t125 = qJD(1) * t96;
	t88 = qJ(2) + pkin(9);
	t86 = sin(t88);
	t115 = t86 * t125;
	t81 = t86 ^ 2;
	t87 = cos(t88);
	t84 = 0.1e1 / t87 ^ 2;
	t132 = t81 * t84;
	t89 = t95 ^ 2;
	t76 = t89 * t132 + 0.1e1;
	t74 = 0.1e1 / t76;
	t83 = 0.1e1 / t87;
	t60 = (-(-t87 * t123 - t115) * t83 + t123 * t132) * t74;
	t144 = t60 - t123;
	t94 = t96 ^ 2;
	t143 = qJD(1) * (0.1e1 / t89 * t94 + 0.1e1) * t129;
	t127 = t95 * t86;
	t73 = atan2(-t127, -t87);
	t71 = sin(t73);
	t117 = t71 * t127;
	t72 = cos(t73);
	t67 = -t72 * t87 - t117;
	t64 = 0.1e1 / t67;
	t65 = 0.1e1 / t67 ^ 2;
	t142 = -0.2e1 * t86;
	t141 = t74 - 0.1e1;
	t134 = t72 * t86;
	t56 = (-t60 * t95 + qJD(2)) * t134 + (t144 * t87 - t115) * t71;
	t140 = t56 * t64 * t65;
	t139 = t60 * t86;
	t138 = t65 * t86;
	t137 = t65 * t96;
	t130 = t83 * t86;
	t80 = t86 * t81;
	t82 = t87 ^ 2;
	t103 = qJD(2) * (t80 * t83 / t82 + t130);
	t108 = t81 * t95 * t125;
	t136 = (t89 * t103 + t84 * t108) / t76 ^ 2;
	t135 = t71 * t95;
	t133 = t81 * t83;
	t131 = t81 * t94;
	t91 = 0.1e1 / t95 ^ 2;
	t128 = t91 * t94;
	t126 = qJD(1) * t95;
	t124 = qJD(2) * t87;
	t122 = qJD(2) * t96;
	t107 = t86 * t94 * t124;
	t63 = t65 * t131 + 0.1e1;
	t121 = 0.2e1 * (-t131 * t140 + (t107 - t108) * t65) / t63 ^ 2;
	t120 = 0.2e1 * t140;
	t79 = t82 * t128 + 0.1e1;
	t119 = 0.2e1 * (-t91 * t107 - t82 * t143) / t79 ^ 2;
	t118 = t86 * t137;
	t116 = t74 * t133;
	t114 = 0.1e1 + t132;
	t113 = 0.1e1 + t128;
	t112 = t86 * t121;
	t111 = t136 * t142;
	t110 = t136 * t145;
	t109 = t95 * t116;
	t106 = t114 * t96;
	t105 = t113 * t86;
	t77 = 0.1e1 / t79;
	t69 = t114 * t95 * t74;
	t61 = 0.1e1 / t63;
	t59 = (t141 * t86 * t71 - t72 * t109) * t96;
	t58 = -t87 * t135 + t134 + (-t72 * t127 + t71 * t87) * t69;
	t57 = -t114 * t110 + (qJD(1) * t106 + t103 * t145) * t74;
	t1 = [t96 * t83 * t111 + (qJD(2) * t106 - t126 * t130) * t74, t57, 0, 0, 0, 0; (t64 * t112 + (-t64 * t124 + (qJD(1) * t59 + t56) * t138) * t61) * t95 + (t65 * t112 * t59 + (-((t60 * t109 + t141 * t124 + t111) * t71 + (t110 * t133 - t139 + (t139 + (-t80 * t84 + t142) * t123) * t74) * t72) * t118 + (t86 * t120 - t65 * t124) * t59 + (-t64 + ((-t89 + t94) * t72 * t116 + t141 * t117) * t65) * t86 * qJD(1)) * t61) * t96, (t58 * t138 - t64 * t87) * t96 * t121 + ((-t64 * t126 + (-qJD(2) * t58 - t56) * t137) * t87 + (-t64 * t122 - (-t57 * t72 * t95 - t144 * t71 + (-qJD(2) * t71 - t125 * t72 + t135 * t60) * t69) * t118 + (t96 * t120 + t65 * t126) * t58 - ((t57 - t125) * t71 + ((-t69 * t95 + 0.1e1) * qJD(2) + (t69 - t95) * t60) * t72) * t87 * t137) * t86) * t61, 0, 0, 0, 0; t113 * t87 * t119 + (qJD(2) * t105 + 0.2e1 * t87 * t143) * t77, t86 * t119 * t129 + (-t87 * t90 * t122 + qJD(1) * t105) * t77, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:54
	% EndTime: 2019-10-10 09:19:55
	% DurationCPUTime: 1.00s
	% Computational Cost: add. (1897->83), mult. (2191->183), div. (456->12), fcn. (2616->9), ass. (0->85)
	t107 = qJ(2) + pkin(9);
	t106 = cos(t107);
	t104 = t106 ^ 2;
	t105 = sin(t107);
	t153 = 0.1e1 / t105 ^ 2 * t104;
	t112 = sin(qJ(1));
	t127 = 0.1e1 + t153;
	t108 = t112 ^ 2;
	t99 = t108 * t153 + 0.1e1;
	t97 = 0.1e1 / t99;
	t124 = t127 * t97;
	t80 = t112 * t124;
	t169 = t112 * t80 - 0.1e1;
	t100 = 0.1e1 / t105;
	t166 = t106 * t153;
	t122 = qJD(2) * (-t106 - t166) * t100;
	t113 = cos(qJ(1));
	t145 = qJD(1) * t113;
	t133 = t112 * t145;
	t168 = (t108 * t122 + t133 * t153) / t99 ^ 2;
	t110 = sin(pkin(10));
	t142 = qJD(2) * t113;
	t131 = t106 * t142;
	t149 = t112 * t110;
	t111 = cos(pkin(10));
	t151 = t111 * t113;
	t93 = -t105 * t149 + t151;
	t85 = t93 * qJD(1) + t110 * t131;
	t148 = t112 * t111;
	t152 = t110 * t113;
	t91 = t105 * t152 + t148;
	t88 = 0.1e1 / t91 ^ 2;
	t167 = t85 * t88;
	t150 = t112 * t106;
	t96 = atan2(-t150, t105);
	t94 = sin(t96);
	t136 = t94 * t150;
	t95 = cos(t96);
	t78 = t105 * t95 - t136;
	t75 = 0.1e1 / t78;
	t87 = 0.1e1 / t91;
	t76 = 0.1e1 / t78 ^ 2;
	t165 = t97 - 0.1e1;
	t143 = qJD(2) * t112;
	t158 = t105 * t94;
	t71 = ((t105 * t143 - t106 * t145) * t100 + t143 * t153) * t97;
	t66 = (-t71 + t143) * t158 + (-t94 * t145 + (-t112 * t71 + qJD(2)) * t95) * t106;
	t164 = t66 * t75 * t76;
	t109 = t113 ^ 2;
	t74 = t104 * t109 * t76 + 0.1e1;
	t72 = 0.1e1 / t74;
	t162 = t72 * t76;
	t161 = t75 * t72;
	t160 = t87 * t167;
	t123 = t105 * t151 - t149;
	t159 = t88 * t123;
	t156 = t113 * t76;
	t155 = qJD(2) * t80;
	t154 = t100 * t104;
	t147 = qJD(1) * t106;
	t146 = qJD(1) * t112;
	t144 = qJD(2) * t105;
	t134 = t76 * t144;
	t141 = 0.2e1 * (-t109 * t106 * t134 + (-t109 * t164 - t76 * t133) * t104) / t74 ^ 2;
	t140 = 0.2e1 * t164;
	t86 = t123 ^ 2;
	t83 = t86 * t88 + 0.1e1;
	t92 = t105 * t148 + t152;
	t84 = t92 * qJD(1) - t111 * t131;
	t139 = 0.2e1 * (-t84 * t159 - t86 * t160) / t83 ^ 2;
	t138 = 0.2e1 * t168;
	t137 = t100 * t112 * t97;
	t135 = t106 * t165;
	t132 = t106 * t143;
	t130 = t75 * t141;
	t129 = t76 * t141;
	t128 = -0.2e1 * t123 * t160;
	t126 = t106 * t138;
	t125 = t104 * t137;
	t81 = 0.1e1 / t83;
	t121 = (-t110 * t159 + t111 * t87) * t81;
	t70 = (t95 * t125 + t94 * t135) * t113;
	t68 = -t169 * t95 * t106 + (t112 - t80) * t158;
	t67 = t124 * t145 + 0.2e1 * (t122 * t97 - t127 * t168) * t112;
	t1 = [t137 * t147 + (qJD(2) * t124 + t100 * t126) * t113, t67, 0, 0, 0, 0; (t144 * t161 + (t130 + (qJD(1) * t70 + t66) * t162) * t106) * t112 + (t70 * t129 * t106 + (t70 * t134 + (t70 * t140 + ((t71 * t125 + t165 * t144 + t126) * t94 + (-t71 * t135 + (t138 * t154 + (0.2e1 * t106 + t166) * t97 * qJD(2)) * t112) * t95) * t156) * t106 + (-t75 + (t165 * t136 - (-t108 + t109) * t97 * t95 * t154) * t76) * t147) * t72) * t113, (t146 * t161 + (t130 + (qJD(2) * t68 + t66) * t162) * t113) * t105 + (t68 * t113 * t129 + (-t75 * t142 + (t113 * t140 + t76 * t146) * t68 + (-((-t112 * t67 - t145 * t80) * t95 + (t169 * t71 + t143 - t155) * t94) * t106 - ((-t67 + t145) * t94 + (-t71 * t80 - qJD(2) + (t71 + t155) * t112) * t95) * t105) * t156) * t72) * t106, 0, 0, 0, 0; (-t93 * t159 - t87 * t92) * t139 + ((t123 * qJD(1) + t111 * t132) * t87 + t93 * t128 + (-t92 * t85 + (-t91 * qJD(1) - t110 * t132) * t123 - t93 * t84) * t88) * t81, t105 * t121 * t142 + (t121 * t146 + ((t87 * t139 + t81 * t167) * t111 + (-t139 * t159 + (-t84 * t88 + t128) * t81) * t110) * t113) * t106, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:54
	% EndTime: 2019-10-10 09:19:55
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (2429->94), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->95)
	t146 = qJ(2) + pkin(9);
	t142 = sin(t146);
	t137 = 0.1e1 / t142 ^ 2;
	t144 = cos(t146);
	t140 = t144 ^ 2;
	t194 = t137 * t140;
	t213 = t144 * t194;
	t149 = sin(qJ(1));
	t170 = 0.1e1 + t194;
	t212 = t149 * t170;
	t166 = qJD(1) * t142 + qJD(6);
	t150 = cos(qJ(1));
	t182 = qJD(2) * t150;
	t211 = -t144 * t182 + t166 * t149;
	t183 = qJD(2) * t149;
	t210 = t144 * t183 + t166 * t150;
	t187 = t149 * t144;
	t131 = atan2(-t187, t142);
	t130 = cos(t131);
	t129 = sin(t131);
	t175 = t129 * t187;
	t115 = t130 * t142 - t175;
	t112 = 0.1e1 / t115;
	t145 = pkin(10) + qJ(6);
	t141 = sin(t145);
	t143 = cos(t145);
	t188 = t149 * t143;
	t190 = t142 * t150;
	t126 = t141 * t190 + t188;
	t122 = 0.1e1 / t126;
	t136 = 0.1e1 / t142;
	t113 = 0.1e1 / t115 ^ 2;
	t123 = 0.1e1 / t126 ^ 2;
	t147 = t149 ^ 2;
	t134 = t147 * t194 + 0.1e1;
	t132 = 0.1e1 / t134;
	t209 = t132 - 0.1e1;
	t185 = qJD(1) * t150;
	t173 = t144 * t185;
	t106 = ((t142 * t183 - t173) * t136 + t183 * t194) * t132;
	t196 = t130 * t144;
	t101 = (-t106 * t149 + qJD(2)) * t196 + (-t173 + (-t106 + t183) * t142) * t129;
	t208 = t101 * t112 * t113;
	t167 = qJD(6) * t142 + qJD(1);
	t161 = t167 * t150;
	t107 = t141 * t161 + t211 * t143;
	t189 = t143 * t150;
	t125 = t141 * t149 - t142 * t189;
	t121 = t125 ^ 2;
	t120 = t121 * t123 + 0.1e1;
	t198 = t123 * t125;
	t108 = -t211 * t141 + t143 * t161;
	t204 = t108 * t122 * t123;
	t207 = (t107 * t198 - t121 * t204) / t120 ^ 2;
	t206 = t106 * t129;
	t205 = t106 * t144;
	t203 = t113 * t144;
	t202 = t113 * t150;
	t195 = t136 * t144;
	t159 = qJD(2) * (-t136 * t213 - t195);
	t192 = t140 * t149;
	t164 = t185 * t192;
	t201 = (t137 * t164 + t147 * t159) / t134 ^ 2;
	t119 = t132 * t212;
	t200 = t119 * t149;
	t199 = t122 * t143;
	t197 = t125 * t141;
	t148 = t150 ^ 2;
	t193 = t140 * t148;
	t191 = t142 * t149;
	t186 = qJD(1) * t149;
	t184 = qJD(2) * t142;
	t111 = t113 * t193 + 0.1e1;
	t181 = 0.2e1 * (-t193 * t208 + (-t144 * t148 * t184 - t164) * t113) / t111 ^ 2;
	t180 = 0.2e1 * t208;
	t179 = 0.2e1 * t207;
	t178 = -0.2e1 * t201;
	t177 = t144 * t202;
	t176 = t144 * t201;
	t174 = t136 * t192;
	t169 = t144 * t181;
	t168 = 0.2e1 * t125 * t204;
	t165 = t132 * t174;
	t163 = t170 * t150;
	t162 = t167 * t149;
	t160 = t123 * t197 + t199;
	t158 = t160 * t150;
	t128 = -t141 * t191 + t189;
	t127 = t141 * t150 + t142 * t188;
	t117 = 0.1e1 / t120;
	t109 = 0.1e1 / t111;
	t105 = (t209 * t144 * t129 + t130 * t165) * t150;
	t104 = t129 * t191 + t196 + (-t129 * t142 - t130 * t187) * t119;
	t102 = t178 * t212 + (qJD(1) * t163 + 0.2e1 * t149 * t159) * t132;
	t1 = [0.2e1 * t136 * t150 * t176 + (qJD(2) * t163 + t186 * t195) * t132, t102, 0, 0, 0, 0; (t112 * t169 + (t112 * t184 + (qJD(1) * t105 + t101) * t203) * t109) * t149 + (t113 * t169 * t105 + (-((-t106 * t165 - t209 * t184 - 0.2e1 * t176) * t129 + (t174 * t178 - t205 + (t205 + (-0.2e1 * t144 - t213) * t183) * t132) * t130) * t177 + (t113 * t184 + t144 * t180) * t105 + (-t112 + ((t147 - t148) * t140 * t136 * t132 * t130 + t209 * t175) * t113) * t144 * qJD(1)) * t109) * t150, (t104 * t203 + t112 * t142) * t150 * t181 + ((t112 * t186 + (qJD(2) * t104 + t101) * t202) * t142 + (-t112 * t182 - (-t102 * t130 * t149 + t129 * t183 + t200 * t206 - t206 + (-qJD(2) * t129 - t130 * t185) * t119) * t177 + (t113 * t186 + t150 * t180) * t104 - ((-t102 + t185) * t129 + ((-0.1e1 + t200) * qJD(2) + (-t119 + t149) * t106) * t130) * t113 * t190) * t144) * t109, 0, 0, 0, 0; (-t122 * t127 + t128 * t198) * t179 + (t128 * t168 - t122 * t141 * t162 + t210 * t199 + (t125 * t143 * t162 - t128 * t107 - t127 * t108 + t210 * t197) * t123) * t117, t144 * t158 * t179 + (t158 * t184 + (t160 * t186 + ((qJD(6) * t122 + t168) * t141 + (-t107 * t141 + (-qJD(6) * t125 + t108) * t143) * t123) * t150) * t144) * t117, 0, 0, 0, -0.2e1 * t207 + 0.2e1 * (t107 * t117 * t123 + (-t117 * t204 - t123 * t207) * t125) * t125;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end