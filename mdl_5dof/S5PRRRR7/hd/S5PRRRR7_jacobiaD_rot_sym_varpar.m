% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRR7
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
%   Wie in S5PRRRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:37
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRRR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:43
	% EndTime: 2019-10-24 10:37:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:43
	% EndTime: 2019-10-24 10:37:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:42
	% EndTime: 2019-10-24 10:37:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:43
	% EndTime: 2019-10-24 10:37:43
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (429->41), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->59)
	t97 = sin(qJ(2));
	t90 = t97 ^ 2;
	t99 = cos(qJ(2));
	t123 = t90 / t99 ^ 2;
	t111 = 0.1e1 + t123;
	t94 = sin(pkin(9));
	t135 = t111 * t94;
	t134 = qJD(2) * t97;
	t95 = cos(pkin(9));
	t120 = t95 * t99;
	t96 = sin(qJ(3));
	t98 = cos(qJ(3));
	t77 = t98 * t120 + t94 * t96;
	t121 = t94 * t97;
	t80 = atan2(-t121, -t99);
	t78 = sin(t80);
	t79 = cos(t80);
	t63 = -t78 * t121 - t79 * t99;
	t60 = 0.1e1 / t63;
	t73 = 0.1e1 / t77;
	t124 = t78 * t99;
	t108 = -t94 * t124 + t79 * t97;
	t115 = t79 * t121;
	t109 = -t115 + t124;
	t87 = t94 ^ 2;
	t83 = t87 * t123 + 0.1e1;
	t81 = 0.1e1 / t83;
	t67 = t81 * t135;
	t54 = t109 * t67 + t108;
	t132 = 0.2e1 * t54;
	t61 = 0.1e1 / t63 ^ 2;
	t74 = 0.1e1 / t77 ^ 2;
	t118 = qJD(2) * t99;
	t127 = t61 * t97;
	t64 = qJD(2) * t67;
	t53 = t108 * qJD(2) + t109 * t64;
	t130 = t53 * t60 * t61;
	t88 = t95 ^ 2;
	t59 = t88 * t90 * t61 + 0.1e1;
	t131 = (t118 * t127 - t90 * t130) * t88 / t59 ^ 2;
	t76 = t96 * t120 - t94 * t98;
	t125 = t74 * t76;
	t113 = t95 * t134;
	t117 = qJD(3) * t76;
	t70 = -t98 * t113 - t117;
	t126 = t70 * t73 * t74;
	t72 = t76 ^ 2;
	t68 = t72 * t74 + 0.1e1;
	t69 = -t77 * qJD(3) + t96 * t113;
	t129 = (-t69 * t125 - t72 * t126) / t68 ^ 2;
	t57 = 0.1e1 / t59;
	t128 = t57 * t61;
	t119 = t67 - t94;
	t116 = -0.2e1 * t129;
	t112 = t67 * t94 - 0.1e1;
	t110 = t98 * t125 - t73 * t96;
	t65 = 0.1e1 / t68;
	t56 = 0.2e1 * (t81 - t111 / t83 ^ 2 * t87) / t99 * t134 * t135;
	t1 = [0, t56, 0, 0, 0; 0, ((-0.2e1 * t60 * t131 + (-qJD(2) * t54 - t53) * t128) * t99 + (t61 * t131 * t132 + (t56 * t61 * t115 + t130 * t132 - qJD(2) * t60 - (-t119 * qJD(2) + t112 * t64) * t78 * t127) * t57 - (t56 * t78 + (-t112 * qJD(2) + t119 * t64) * t79) * t99 * t128) * t97) * t95, 0, 0, 0; 0, (t110 * t97 * t116 + (t110 * t118 + ((-qJD(3) * t73 - 0.2e1 * t76 * t126) * t98 + (-t69 * t98 + (t70 - t117) * t96) * t74) * t97) * t65) * t95, t116 + 0.2e1 * (-t65 * t69 * t74 + (-t65 * t126 - t74 * t129) * t76) * t76, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:43
	% EndTime: 2019-10-24 10:37:43
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (788->45), mult. (1346->105), div. (247->12), fcn. (1511->9), ass. (0->61)
	t125 = sin(pkin(9));
	t127 = sin(qJ(2));
	t120 = t127 ^ 2;
	t128 = cos(qJ(2));
	t149 = t120 / t128 ^ 2;
	t138 = 0.1e1 + t149;
	t163 = t125 * t138;
	t162 = qJD(2) * t127;
	t124 = qJ(3) + qJ(4);
	t114 = sin(t124);
	t115 = cos(t124);
	t126 = cos(pkin(9));
	t146 = t126 * t128;
	t104 = t125 * t114 + t115 * t146;
	t147 = t125 * t127;
	t108 = atan2(-t147, -t128);
	t106 = sin(t108);
	t107 = cos(t108);
	t95 = -t106 * t147 - t107 * t128;
	t92 = 0.1e1 / t95;
	t116 = t125 ^ 2;
	t112 = t116 * t149 + 0.1e1;
	t110 = 0.1e1 / t112;
	t97 = t110 * t163;
	t140 = t125 * t97 - 0.1e1;
	t150 = t107 * t127;
	t151 = t106 * t128;
	t153 = t125 - t97;
	t82 = -t140 * t150 - t153 * t151;
	t160 = 0.2e1 * t82;
	t100 = 0.1e1 / t104;
	t101 = 0.1e1 / t104 ^ 2;
	t93 = 0.1e1 / t95 ^ 2;
	t117 = t126 ^ 2;
	t145 = qJD(2) * t128;
	t154 = t127 * t93;
	t143 = t107 * t147;
	t96 = qJD(2) * t97;
	t81 = (-t143 + t151) * t96 + (-t125 * t151 + t150) * qJD(2);
	t157 = t81 * t92 * t93;
	t86 = t117 * t120 * t93 + 0.1e1;
	t159 = (-t120 * t157 + t145 * t154) * t117 / t86 ^ 2;
	t142 = t114 * t146;
	t103 = -t125 * t115 + t142;
	t152 = t101 * t103;
	t118 = qJD(3) + qJD(4);
	t139 = t126 * t162;
	t91 = -t118 * t142 + (t118 * t125 - t139) * t115;
	t155 = t100 * t101 * t91;
	t99 = t103 ^ 2;
	t89 = t99 * t101 + 0.1e1;
	t90 = -t104 * t118 + t114 * t139;
	t158 = (-t90 * t152 - t99 * t155) / t89 ^ 2;
	t84 = 0.1e1 / t86;
	t156 = t84 * t93;
	t144 = -0.2e1 * t158;
	t137 = -t100 * t114 + t115 * t152;
	t87 = 0.1e1 / t89;
	t83 = 0.2e1 * (t110 - t138 * t116 / t112 ^ 2) / t128 * t162 * t163;
	t78 = t144 + 0.2e1 * (-t101 * t87 * t90 + (-t101 * t158 - t87 * t155) * t103) * t103;
	t1 = [0, t83, 0, 0, 0; 0, ((-0.2e1 * t92 * t159 + (-qJD(2) * t82 - t81) * t156) * t128 + (t93 * t159 * t160 + (t83 * t93 * t143 + t157 * t160 - qJD(2) * t92 - (t153 * qJD(2) + t140 * t96) * t106 * t154) * t84 - (t106 * t83 + (qJD(2) + (-0.2e1 * t125 + t97) * t96) * t107) * t128 * t156) * t127) * t126, 0, 0, 0; 0, (t137 * t87 * t145 + (t137 * t144 + ((-t100 * t118 - 0.2e1 * t103 * t155) * t115 + (-t115 * t90 + (-t103 * t118 + t91) * t114) * t101) * t87) * t127) * t126, t78, t78, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:43
	% EndTime: 2019-10-24 10:37:43
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (1329->44), mult. (1526->104), div. (265->12), fcn. (1707->9), ass. (0->62)
	t130 = sin(pkin(9));
	t132 = sin(qJ(2));
	t126 = t132 ^ 2;
	t133 = cos(qJ(2));
	t157 = t126 / t133 ^ 2;
	t145 = 0.1e1 + t157;
	t169 = t130 * t145;
	t168 = qJD(2) * t132;
	t122 = qJ(3) + qJ(4) + qJ(5);
	t119 = sin(t122);
	t120 = cos(t122);
	t131 = cos(pkin(9));
	t154 = t131 * t133;
	t108 = t130 * t119 + t120 * t154;
	t123 = t130 ^ 2;
	t117 = t123 * t157 + 0.1e1;
	t115 = 0.1e1 / t117;
	t102 = t115 * t169;
	t155 = t130 * t132;
	t114 = atan2(-t155, -t133);
	t113 = cos(t114);
	t112 = sin(t114);
	t158 = t112 * t133;
	t142 = t113 * t132 - t130 * t158;
	t150 = t113 * t155;
	t143 = -t150 + t158;
	t87 = t143 * t102 + t142;
	t166 = 0.2e1 * t87;
	t100 = -t112 * t155 - t113 * t133;
	t97 = 0.1e1 / t100;
	t104 = 0.1e1 / t108;
	t98 = 0.1e1 / t100 ^ 2;
	t105 = 0.1e1 / t108 ^ 2;
	t124 = t131 ^ 2;
	t152 = qJD(2) * t133;
	t160 = t132 * t98;
	t101 = qJD(2) * t102;
	t86 = t142 * qJD(2) + t143 * t101;
	t163 = t86 * t97 * t98;
	t96 = t124 * t126 * t98 + 0.1e1;
	t165 = (-t126 * t163 + t152 * t160) * t124 / t96 ^ 2;
	t149 = t119 * t154;
	t107 = -t130 * t120 + t149;
	t103 = t107 ^ 2;
	t159 = t105 * t107;
	t121 = qJD(3) + qJD(4) + qJD(5);
	t147 = t131 * t168;
	t95 = -t121 * t149 + (t121 * t130 - t147) * t120;
	t161 = t104 * t105 * t95;
	t91 = t103 * t105 + 0.1e1;
	t94 = -t108 * t121 + t119 * t147;
	t164 = (-t103 * t161 - t94 * t159) / t91 ^ 2;
	t92 = 0.1e1 / t96;
	t162 = t92 * t98;
	t153 = t102 - t130;
	t151 = -0.2e1 * t164;
	t146 = t102 * t130 - 0.1e1;
	t144 = -t104 * t119 + t120 * t159;
	t89 = 0.1e1 / t91;
	t88 = 0.2e1 * (t115 - t145 / t117 ^ 2 * t123) / t133 * t168 * t169;
	t83 = t151 + 0.2e1 * (-t105 * t89 * t94 + (-t105 * t164 - t89 * t161) * t107) * t107;
	t1 = [0, t88, 0, 0, 0; 0, ((-0.2e1 * t97 * t165 + (-qJD(2) * t87 - t86) * t162) * t133 + (t98 * t165 * t166 + (t88 * t98 * t150 + t163 * t166 - qJD(2) * t97 - (-t153 * qJD(2) + t146 * t101) * t112 * t160) * t92 - (t112 * t88 + (-t146 * qJD(2) + t153 * t101) * t113) * t133 * t162) * t132) * t131, 0, 0, 0; 0, (t144 * t89 * t152 + (t144 * t151 + ((-t104 * t121 - 0.2e1 * t107 * t161) * t120 + (-t120 * t94 + (-t107 * t121 + t95) * t119) * t105) * t89) * t132) * t131, t83, t83, t83;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end