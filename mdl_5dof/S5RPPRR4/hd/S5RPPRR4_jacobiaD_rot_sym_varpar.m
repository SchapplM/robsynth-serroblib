% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPRR4
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
%   Wie in S5RPPRR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPPRR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
JaD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:42
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (265->32), mult. (614->97), div. (108->12), fcn. (792->9), ass. (0->49)
	t86 = sin(pkin(8));
	t79 = t86 ^ 2;
	t88 = cos(pkin(8));
	t81 = 0.1e1 / t88 ^ 2;
	t89 = sin(qJ(1));
	t83 = t89 ^ 2;
	t77 = t79 * t81 * t83 + 0.1e1;
	t90 = cos(qJ(1));
	t84 = t90 ^ 2;
	t106 = 0.1e1 / t77 ^ 2 * t84;
	t112 = t106 * t81;
	t75 = 0.1e1 / t77;
	t111 = (t75 - 0.1e1) * t86;
	t104 = t86 * t89;
	t74 = atan2(-t104, -t88);
	t72 = sin(t74);
	t73 = cos(t74);
	t58 = -t72 * t104 - t73 * t88;
	t55 = 0.1e1 / t58;
	t85 = sin(pkin(9));
	t102 = t89 * t85;
	t103 = t88 * t90;
	t87 = cos(pkin(9));
	t71 = t87 * t103 + t102;
	t65 = 0.1e1 / t71;
	t80 = 0.1e1 / t88;
	t56 = 0.1e1 / t58 ^ 2;
	t66 = 0.1e1 / t71 ^ 2;
	t109 = t56 * t90;
	t101 = t89 * t87;
	t70 = t85 * t103 - t101;
	t108 = t66 * t70;
	t69 = -t88 * t101 + t85 * t90;
	t107 = t69 * t70;
	t105 = t79 * t80;
	t100 = qJD(1) * t89;
	t99 = t80 * t112;
	t68 = -t88 * t102 - t87 * t90;
	t51 = (-t73 * t75 * t89 * t105 + t72 * t111) * t90;
	t78 = t86 * t79;
	t67 = t65 * t66;
	t64 = t70 ^ 2;
	t63 = t69 * qJD(1);
	t62 = t68 * qJD(1);
	t61 = t64 * t66 + 0.1e1;
	t57 = t55 * t56;
	t54 = t56 * t79 * t84 + 0.1e1;
	t50 = qJD(1) * t51;
	t1 = [(-t75 * t80 * t86 - 0.2e1 * t78 * t99) * t100, 0, 0, 0, 0; (0.2e1 * (t51 * t109 + t55 * t89) / t54 ^ 2 * (-t50 * t57 * t84 - t100 * t109) * t79 + ((0.2e1 * t51 * t57 * t90 + t56 * t89) * t50 + (-t90 * t55 + ((t51 + (t78 * t112 + t111) * t90 * t72) * t89 - (0.2e1 * t83 * t79 ^ 2 * t99 + (t106 + (t83 - 0.2e1 * t84) * t75) * t105) * t90 * t73) * t56) * qJD(1)) / t54) * t86, 0, 0, 0, 0; 0.2e1 * (t66 * t107 - t65 * t68) / t61 ^ 2 * (-t63 * t64 * t67 + t62 * t108) + (-t69 * t62 * t66 + (0.2e1 * t67 * t107 - t68 * t66) * t63 + (t71 * t108 - t70 * t65) * qJD(1)) / t61, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:42
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (537->40), mult. (853->110), div. (126->12), fcn. (1047->9), ass. (0->55)
	t107 = cos(pkin(8));
	t101 = 0.1e1 / t107 ^ 2;
	t108 = sin(qJ(1));
	t104 = t108 ^ 2;
	t106 = sin(pkin(8));
	t99 = t106 ^ 2;
	t93 = t99 * t104 * t101 + 0.1e1;
	t92 = 0.1e1 / t93 ^ 2;
	t137 = t101 * t92;
	t109 = cos(qJ(1));
	t123 = t107 * t109;
	t103 = pkin(9) + qJ(4);
	t96 = sin(t103);
	t97 = cos(t103);
	t86 = t108 * t96 + t97 * t123;
	t125 = t106 * t108;
	t90 = atan2(-t125, -t107);
	t87 = sin(t90);
	t120 = t87 * t125;
	t88 = cos(t90);
	t78 = -t88 * t107 - t120;
	t75 = 0.1e1 / t78;
	t80 = 0.1e1 / t86;
	t100 = 0.1e1 / t107;
	t76 = 0.1e1 / t78 ^ 2;
	t81 = 0.1e1 / t86 ^ 2;
	t91 = 0.1e1 / t93;
	t136 = t91 - 0.1e1;
	t124 = t107 * t108;
	t84 = t109 * t96 - t97 * t124;
	t85 = -t108 * t97 + t96 * t123;
	t71 = t84 * qJD(1) - t85 * qJD(4);
	t133 = t71 * t80 * t81;
	t118 = t109 * t97 + t96 * t124;
	t70 = t118 * qJD(1) - t86 * qJD(4);
	t134 = t70 * t81;
	t79 = t85 ^ 2;
	t74 = t79 * t81 + 0.1e1;
	t135 = (-t79 * t133 - t85 * t134) / t74 ^ 2;
	t132 = t84 * t85;
	t131 = t100 * t99;
	t130 = t100 * t137;
	t105 = t109 ^ 2;
	t129 = t105 * t76;
	t128 = t108 * t76;
	t126 = t109 * t76;
	t122 = qJD(1) * t108;
	t119 = t88 * t91 * t131;
	t66 = (t136 * t87 * t106 - t108 * t119) * t109;
	t98 = t106 * t99;
	t77 = t75 * t76;
	t72 = 0.1e1 / t74;
	t69 = t99 * t129 + 0.1e1;
	t65 = qJD(1) * t66;
	t1 = [(-t100 * t106 * t91 - 0.2e1 * t105 * t98 * t130) * t122, 0, 0, 0, 0; (0.2e1 * (t108 * t75 + t66 * t126) / t69 ^ 2 * (-t105 * t65 * t77 - t122 * t126) * t99 + ((0.2e1 * t109 * t66 * t77 + t128) * t65 + (t66 * t128 + (-t75 - (-t108 * t87 * t98 * t137 + (0.2e1 * t104 * t99 ^ 2 * t130 + (-0.2e1 * t91 + t92) * t131) * t88) * t129 + (-t104 * t119 + t136 * t120) * t76) * t109) * qJD(1)) / t69) * t106, 0, 0, 0, 0; 0.2e1 * (t118 * t80 + t81 * t132) * t135 + ((-t85 * qJD(1) + t84 * qJD(4)) * t80 + 0.2e1 * t132 * t133 + (t118 * t71 - (-t86 * qJD(1) + t118 * qJD(4)) * t85 + t84 * t70) * t81) * t72, 0, 0, -0.2e1 * t135 + 0.2e1 * (-t72 * t134 + (-t72 * t133 - t81 * t135) * t85) * t85, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:42
	% EndTime: 2022-01-23 09:17:42
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (1005->41), mult. (1047->107), div. (144->12), fcn. (1257->9), ass. (0->55)
	t134 = cos(pkin(8));
	t128 = 0.1e1 / t134 ^ 2;
	t133 = sin(pkin(8));
	t126 = t133 ^ 2;
	t135 = sin(qJ(1));
	t131 = t135 ^ 2;
	t120 = t131 * t126 * t128 + 0.1e1;
	t136 = cos(qJ(1));
	t132 = t136 ^ 2;
	t155 = 0.1e1 / t120 ^ 2 * t132;
	t163 = t128 * t155;
	t118 = 0.1e1 / t120;
	t162 = (t118 - 0.1e1) * t133;
	t124 = pkin(9) + qJ(4) + qJ(5);
	t123 = cos(t124);
	t149 = t136 * t123;
	t122 = sin(t124);
	t153 = t135 * t122;
	t112 = t134 * t149 + t153;
	t151 = t135 * t133;
	t117 = atan2(-t151, -t134);
	t114 = sin(t117);
	t115 = cos(t117);
	t104 = -t114 * t151 - t115 * t134;
	t101 = 0.1e1 / t104;
	t106 = 0.1e1 / t112;
	t127 = 0.1e1 / t134;
	t102 = 0.1e1 / t104 ^ 2;
	t107 = 0.1e1 / t112 ^ 2;
	t150 = t136 * t122;
	t152 = t135 * t123;
	t111 = t134 * t150 - t152;
	t105 = t111 ^ 2;
	t130 = qJD(4) + qJD(5);
	t145 = t134 * t153 + t149;
	t93 = t145 * qJD(1) - t112 * t130;
	t158 = t93 * t107;
	t110 = -t134 * t152 + t150;
	t94 = t110 * qJD(1) - t111 * t130;
	t159 = t106 * t107 * t94;
	t97 = t105 * t107 + 0.1e1;
	t160 = (-t105 * t159 - t111 * t158) / t97 ^ 2;
	t157 = t102 * t136;
	t156 = t110 * t111;
	t154 = t126 * t127;
	t148 = qJD(1) * t135;
	t147 = t127 * t163;
	t92 = (-t115 * t118 * t135 * t154 + t114 * t162) * t136;
	t125 = t133 * t126;
	t103 = t101 * t102;
	t100 = t132 * t126 * t102 + 0.1e1;
	t95 = 0.1e1 / t97;
	t91 = qJD(1) * t92;
	t88 = -0.2e1 * t160 + 0.2e1 * (-t95 * t158 + (-t107 * t160 - t95 * t159) * t111) * t111;
	t1 = [(-t118 * t127 * t133 - 0.2e1 * t125 * t147) * t148, 0, 0, 0, 0; (0.2e1 * (t101 * t135 + t92 * t157) / t100 ^ 2 * (-t103 * t132 * t91 - t148 * t157) * t126 + ((0.2e1 * t103 * t136 * t92 + t102 * t135) * t91 + (-t136 * t101 + ((t92 + (t125 * t163 + t162) * t136 * t114) * t135 - (0.2e1 * t131 * t126 ^ 2 * t147 + (t155 + (t131 - 0.2e1 * t132) * t118) * t154) * t136 * t115) * t102) * qJD(1)) / t100) * t133, 0, 0, 0, 0; 0.2e1 * (t106 * t145 + t107 * t156) * t160 + ((-t111 * qJD(1) + t110 * t130) * t106 + 0.2e1 * t156 * t159 + (t145 * t94 - (-t112 * qJD(1) + t145 * t130) * t111 + t110 * t93) * t107) * t95, 0, 0, t88, t88;];
	JaD_rot = t1;
end