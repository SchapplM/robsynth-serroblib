% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR7
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
%   Wie in S6RPRPPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:43
	% DurationCPUTime: 0.68s
	% Computational Cost: add. (1694->69), mult. (1839->158), div. (436->14), fcn. (2165->7), ass. (0->74)
	t96 = cos(qJ(1));
	t144 = 0.2e1 * t96;
	t92 = 0.1e1 / t96;
	t95 = sin(qJ(1));
	t126 = t92 * t95;
	t90 = t95 ^ 2;
	t91 = t96 ^ 2;
	t104 = qJD(1) * (t90 / t91 + 0.1e1) * t126;
	t88 = qJ(3) + pkin(9);
	t86 = sin(t88);
	t122 = qJD(3) * t86;
	t87 = cos(t88);
	t107 = t87 * t90 * t122;
	t93 = 0.1e1 / t96 ^ 2;
	t127 = t90 * t93;
	t80 = t86 ^ 2;
	t79 = t80 * t127 + 0.1e1;
	t143 = -0.2e1 * (t80 * t104 + t93 * t107) / t79 ^ 2;
	t120 = qJD(3) * t96;
	t124 = qJD(1) * t95;
	t114 = t87 * t124;
	t82 = 0.1e1 / t86 ^ 2;
	t85 = t87 ^ 2;
	t129 = t82 * t85;
	t76 = t91 * t129 + 0.1e1;
	t74 = 0.1e1 / t76;
	t81 = 0.1e1 / t86;
	t60 = ((t86 * t120 + t114) * t81 + t120 * t129) * t74;
	t142 = -t60 + t120;
	t125 = t96 * t87;
	t73 = atan2(-t125, t86);
	t71 = sin(t73);
	t72 = cos(t73);
	t67 = -t71 * t125 + t72 * t86;
	t64 = 0.1e1 / t67;
	t65 = 0.1e1 / t67 ^ 2;
	t141 = t74 - 0.1e1;
	t123 = qJD(1) * t96;
	t108 = t85 * t95 * t123;
	t128 = t85 * t90;
	t132 = t72 * t87;
	t56 = (-t60 * t96 + qJD(3)) * t132 + (t142 * t86 + t114) * t71;
	t139 = t56 * t64 * t65;
	t63 = t65 * t128 + 0.1e1;
	t140 = (-t128 * t139 + (-t107 + t108) * t65) / t63 ^ 2;
	t138 = t60 * t87;
	t137 = t65 * t87;
	t136 = t65 * t95;
	t130 = t81 * t87;
	t84 = t87 * t85;
	t103 = qJD(3) * (-t81 / t80 * t84 - t130);
	t135 = (t91 * t103 - t82 * t108) / t76 ^ 2;
	t133 = t71 * t96;
	t131 = t81 * t85;
	t121 = qJD(3) * t95;
	t119 = -0.2e1 * t139;
	t118 = t87 * t140;
	t117 = t87 * t136;
	t116 = t87 * t135;
	t115 = t74 * t131;
	t113 = 0.1e1 + t129;
	t112 = 0.1e1 + t127;
	t111 = t135 * t144;
	t110 = t96 * t115;
	t109 = t141 * t87 * t71;
	t106 = t113 * t95;
	t105 = t112 * t87;
	t77 = 0.1e1 / t79;
	t69 = t113 * t96 * t74;
	t61 = 0.1e1 / t63;
	t59 = (-t72 * t110 - t109) * t95;
	t58 = t86 * t133 + t132 + (-t72 * t125 - t71 * t86) * t69;
	t57 = -t113 * t111 + (-qJD(1) * t106 + t103 * t144) * t74;
	t1 = [-0.2e1 * t95 * t81 * t116 + (-qJD(3) * t106 + t123 * t130) * t74, 0, t57, 0, 0, 0; (0.2e1 * t64 * t118 + (t64 * t122 + (qJD(1) * t59 + t56) * t137) * t61) * t96 + (-0.2e1 * t65 * t118 * t59 + (((t60 * t110 + t141 * t122 + 0.2e1 * t116) * t71 + (t111 * t131 + t138 + (-t138 + (t82 * t84 + 0.2e1 * t87) * t120) * t74) * t72) * t117 + (t87 * t119 - t65 * t122) * t59 + (t64 + ((t90 - t91) * t72 * t115 - t96 * t109) * t65) * t87 * qJD(1)) * t61) * t95, 0, 0.2e1 * (-t58 * t137 - t64 * t86) * t95 * t140 + ((t64 * t123 + (-qJD(3) * t58 - t56) * t136) * t86 + (t64 * t121 + (-t57 * t72 * t96 + t142 * t71 + (-qJD(3) * t71 + t124 * t72 + t133 * t60) * t69) * t117 + (t95 * t119 + t65 * t123) * t58 + ((-t57 - t124) * t71 + ((t69 * t96 - 0.1e1) * qJD(3) + (-t69 + t96) * t60) * t72) * t86 * t136) * t87) * t61, 0, 0, 0; t112 * t86 * t143 + (qJD(3) * t105 + 0.2e1 * t104 * t86) * t77, 0, t87 * t126 * t143 + (-t86 * t92 * t121 + qJD(1) * t105) * t77, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:43
	% DurationCPUTime: 1.02s
	% Computational Cost: add. (1892->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->94)
	t140 = cos(qJ(1));
	t136 = t140 ^ 2;
	t134 = qJ(3) + pkin(9);
	t132 = sin(t134);
	t128 = t132 ^ 2;
	t133 = cos(t134);
	t130 = 0.1e1 / t133 ^ 2;
	t185 = t128 * t130;
	t124 = t136 * t185 + 0.1e1;
	t127 = t132 * t128;
	t129 = 0.1e1 / t133;
	t183 = t129 * t132;
	t148 = qJD(3) * (t127 * t129 * t130 + t183);
	t138 = sin(qJ(1));
	t175 = qJD(1) * t140;
	t152 = t128 * t138 * t175;
	t192 = (-t130 * t152 + t136 * t148) / t124 ^ 2;
	t203 = -0.2e1 * t192;
	t162 = 0.1e1 + t185;
	t202 = t140 * t162;
	t155 = qJD(1) * t133 + qJD(6);
	t172 = qJD(3) * t140;
	t201 = t132 * t172 + t155 * t138;
	t177 = t140 * t132;
	t123 = atan2(t177, t133);
	t119 = sin(t123);
	t120 = cos(t123);
	t105 = t119 * t177 + t120 * t133;
	t102 = 0.1e1 / t105;
	t139 = cos(qJ(6));
	t178 = t139 * t140;
	t137 = sin(qJ(6));
	t180 = t138 * t137;
	t118 = -t133 * t180 + t178;
	t112 = 0.1e1 / t118;
	t103 = 0.1e1 / t105 ^ 2;
	t113 = 0.1e1 / t118 ^ 2;
	t121 = 0.1e1 / t124;
	t200 = t121 - 0.1e1;
	t135 = t138 ^ 2;
	t174 = qJD(3) * t133;
	t184 = t128 * t135;
	t164 = t130 * t172;
	t176 = qJD(1) * t138;
	t165 = t132 * t176;
	t96 = ((t133 * t172 - t165) * t129 + t128 * t164) * t121;
	t157 = -t96 + t172;
	t158 = t140 * t96 - qJD(3);
	t187 = t120 * t132;
	t91 = t158 * t187 + (t157 * t133 - t165) * t119;
	t197 = t102 * t103 * t91;
	t99 = t103 * t184 + 0.1e1;
	t199 = (-t184 * t197 + (t132 * t135 * t174 + t152) * t103) / t99 ^ 2;
	t97 = 0.1e1 / t99;
	t198 = t103 * t97;
	t156 = qJD(6) * t133 + qJD(1);
	t166 = t133 * t178;
	t100 = -qJD(1) * t166 - qJD(6) * t178 + (qJD(3) * t132 * t139 + t156 * t137) * t138;
	t179 = t138 * t139;
	t181 = t137 * t140;
	t117 = t133 * t179 + t181;
	t111 = t117 ^ 2;
	t110 = t111 * t113 + 0.1e1;
	t190 = t113 * t117;
	t173 = qJD(3) * t138;
	t101 = -t156 * t179 + (t132 * t173 - t155 * t140) * t137;
	t194 = t101 * t112 * t113;
	t196 = 0.1e1 / t110 ^ 2 * (-t100 * t190 - t111 * t194);
	t191 = t112 * t139;
	t189 = t117 * t137;
	t188 = t119 * t133;
	t186 = t128 * t129;
	t182 = t132 * t138;
	t171 = 0.2e1 * t199;
	t170 = 0.2e1 * t197;
	t169 = 0.2e1 * t196;
	t168 = t97 * t174;
	t167 = t140 * t186;
	t161 = -0.2e1 * t102 * t199;
	t160 = 0.2e1 * t117 * t194;
	t159 = 0.2e1 * t132 * t192;
	t154 = t121 * t167;
	t153 = t200 * t132 * t119;
	t151 = t162 * t138;
	t150 = t140 * t156;
	t149 = t113 * t189 + t191;
	t116 = -t133 * t181 - t179;
	t115 = t166 - t180;
	t108 = 0.1e1 / t110;
	t107 = t121 * t202;
	t95 = (-t120 * t154 + t153) * t138;
	t93 = t140 * t188 - t187 + (t120 * t177 - t188) * t107;
	t92 = t202 * t203 + (-qJD(1) * t151 + 0.2e1 * t140 * t148) * t121;
	t1 = [t138 * t129 * t159 + (-qJD(3) * t151 - t175 * t183) * t121, 0, t92, 0, 0, 0; (t102 * t168 + (t161 + (-qJD(1) * t95 - t91) * t198) * t132) * t140 + ((-t95 * t168 + (t95 * t171 + ((-t96 * t154 - t200 * t174 + t159) * t119 + (t167 * t203 + t132 * t96 + (t127 * t164 - (t96 - 0.2e1 * t172) * t132) * t121) * t120) * t97 * t138) * t132) * t103 + (t95 * t170 + (-t102 + ((-t135 + t136) * t121 * t120 * t186 - t140 * t153) * t103) * qJD(1)) * t132 * t97) * t138, 0, (t102 * t97 * t175 + (t161 + (-qJD(3) * t93 - t91) * t198) * t138) * t133 + (((-qJD(3) * t102 + t93 * t170) * t138 + (-t93 * t175 + (-(-t107 * t176 + t140 * t92) * t120 - ((-t107 * t140 + 0.1e1) * t96 + (t107 - t140) * qJD(3)) * t119) * t182) * t103) * t97 + (t93 * t171 - ((-t92 - t176) * t119 + (t157 * t107 + t158) * t120) * t97 * t133) * t103 * t138) * t132, 0, 0, 0; (-t112 * t115 + t116 * t190) * t169 + (t116 * t160 - t112 * t137 * t150 - t201 * t191 + (t117 * t139 * t150 + t116 * t100 - t115 * t101 - t201 * t189) * t113) * t108, 0, t149 * t169 * t182 + (-t149 * t133 * t173 + (-t149 * t175 + ((qJD(6) * t112 + t160) * t137 + (t100 * t137 + (-qJD(6) * t117 + t101) * t139) * t113) * t138) * t132) * t108, 0, 0, -0.2e1 * t196 + 0.2e1 * (-t100 * t113 * t108 + (-t108 * t194 - t113 * t196) * t117) * t117;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end