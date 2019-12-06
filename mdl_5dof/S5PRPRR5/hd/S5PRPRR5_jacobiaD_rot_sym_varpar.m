% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPRR5
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
%   Wie in S5PRPRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRPRR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:55:21
	% EndTime: 2019-12-05 15:55:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:55:21
	% EndTime: 2019-12-05 15:55:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:55:21
	% EndTime: 2019-12-05 15:55:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:55:21
	% EndTime: 2019-12-05 15:55:21
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (346->34), mult. (901->88), div. (202->12), fcn. (1017->9), ass. (0->52)
	t80 = sin(pkin(8));
	t83 = sin(qJ(2));
	t75 = t83 ^ 2;
	t84 = cos(qJ(2));
	t102 = t75 / t84 ^ 2;
	t95 = 0.1e1 + t102;
	t113 = t80 * t95;
	t82 = cos(pkin(8));
	t73 = t82 ^ 2;
	t103 = t73 * t75;
	t112 = 0.2e1 * t103;
	t111 = qJD(2) * t83;
	t101 = t80 * t83;
	t68 = atan2(-t101, -t84);
	t66 = sin(t68);
	t67 = cos(t68);
	t53 = -t66 * t101 - t67 * t84;
	t50 = 0.1e1 / t53;
	t100 = t82 * t84;
	t79 = sin(pkin(9));
	t81 = cos(pkin(9));
	t65 = t81 * t100 + t80 * t79;
	t61 = 0.1e1 / t65;
	t51 = 0.1e1 / t53 ^ 2;
	t62 = 0.1e1 / t65 ^ 2;
	t72 = t80 ^ 2;
	t71 = t72 * t102 + 0.1e1;
	t69 = 0.1e1 / t71;
	t58 = t69 * t113;
	t104 = t66 * t84;
	t93 = -t80 * t104 + t67 * t83;
	t98 = t67 * t101;
	t94 = -t98 + t104;
	t44 = t94 * t58 + t93;
	t109 = 0.2e1 * t44;
	t106 = t51 * t84;
	t54 = qJD(2) * t58;
	t43 = t93 * qJD(2) + t94 * t54;
	t107 = t43 * t50 * t51;
	t49 = t51 * t103 + 0.1e1;
	t108 = (t106 * t111 - t75 * t107) * t73 / t49 ^ 2;
	t97 = t79 * t100;
	t64 = -t80 * t81 + t97;
	t105 = t64 * t81;
	t99 = t58 - t80;
	t96 = t58 * t80 - 0.1e1;
	t63 = t61 * t62;
	t60 = t64 ^ 2;
	t57 = t60 * t62 + 0.1e1;
	t47 = 0.1e1 / t49;
	t45 = 0.2e1 * (t69 - t95 / t71 ^ 2 * t72) / t84 * t111 * t113;
	t1 = [0, t45, 0, 0, 0; 0, ((-0.2e1 * t50 * t108 + (-qJD(2) * t44 - t43) * t51 * t47) * t84 + (t51 * t108 * t109 + (t107 * t109 - qJD(2) * t50 - (t45 * t66 + (-t96 * qJD(2) + t99 * t54) * t67) * t106 + (t45 * t98 - (-t99 * qJD(2) + t96 * t54) * t83 * t66) * t51) * t47) * t83) * t82, 0, 0, 0; 0, ((-t62 * t105 + t61 * t79) / t57 ^ 2 * (t60 * t63 * t81 - t62 * t64 * t79) * t112 + (-t61 * t97 + (t63 * t105 * t112 + (t64 * t100 - 0.2e1 * t79 * t103) * t62) * t81) / t57) * qJD(2), 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:55:21
	% EndTime: 2019-12-05 15:55:21
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (614->45), mult. (1166->104), div. (229->12), fcn. (1315->9), ass. (0->57)
	t102 = sin(pkin(8));
	t105 = cos(qJ(2));
	t104 = sin(qJ(2));
	t98 = t104 ^ 2;
	t128 = 0.1e1 / t105 ^ 2 * t98;
	t116 = 0.1e1 + t128;
	t138 = t102 * t116;
	t137 = qJD(2) * t104;
	t103 = cos(pkin(8));
	t121 = t103 * t105;
	t96 = pkin(9) + qJ(4);
	t92 = sin(t96);
	t93 = cos(t96);
	t82 = t102 * t92 + t93 * t121;
	t122 = t102 * t104;
	t86 = atan2(-t122, -t105);
	t83 = sin(t86);
	t84 = cos(t86);
	t71 = -t84 * t105 - t83 * t122;
	t68 = 0.1e1 / t71;
	t78 = 0.1e1 / t82;
	t124 = t84 * t104;
	t125 = t105 * t83;
	t94 = t102 ^ 2;
	t89 = t94 * t128 + 0.1e1;
	t87 = 0.1e1 / t89;
	t75 = t87 * t138;
	t60 = t75 * t125 + t124 + (-t75 * t124 - t125) * t102;
	t135 = 0.2e1 * t60;
	t69 = 0.1e1 / t71 ^ 2;
	t79 = 0.1e1 / t82 ^ 2;
	t120 = qJD(2) * t105;
	t126 = t104 * t69;
	t118 = t84 * t122;
	t72 = qJD(2) * t75;
	t59 = (-t118 + t125) * t72 + (-t102 * t125 + t124) * qJD(2);
	t132 = t59 * t68 * t69;
	t95 = t103 ^ 2;
	t64 = t95 * t98 * t69 + 0.1e1;
	t134 = (t120 * t126 - t98 * t132) * t95 / t64 ^ 2;
	t81 = -t102 * t93 + t92 * t121;
	t129 = t79 * t81;
	t115 = t103 * t137;
	t123 = qJD(4) * t81;
	t74 = -t93 * t115 - t123;
	t130 = t74 * t78 * t79;
	t77 = t81 ^ 2;
	t67 = t77 * t79 + 0.1e1;
	t73 = -t82 * qJD(4) + t92 * t115;
	t133 = (-t73 * t129 - t77 * t130) / t67 ^ 2;
	t62 = 0.1e1 / t64;
	t131 = t62 * t69;
	t119 = -0.2e1 * t133;
	t114 = t93 * t129 - t78 * t92;
	t65 = 0.1e1 / t67;
	t61 = 0.2e1 * (t87 - t116 / t89 ^ 2 * t94) / t105 * t137 * t138;
	t1 = [0, t61, 0, 0, 0; 0, ((-0.2e1 * t68 * t134 + (-qJD(2) * t60 - t59) * t131) * t105 + (t69 * t134 * t135 + (t61 * t69 * t118 + t132 * t135 - qJD(2) * t68 - ((t102 * t75 - 0.1e1) * t72 + (t102 - t75) * qJD(2)) * t83 * t126) * t62 - (t61 * t83 + (qJD(2) + (-0.2e1 * t102 + t75) * t72) * t84) * t105 * t131) * t104) * t103, 0, 0, 0; 0, (t114 * t65 * t120 + (t114 * t119 + ((-qJD(4) * t78 - 0.2e1 * t81 * t130) * t93 + (-t73 * t93 + (t74 - t123) * t92) * t79) * t65) * t104) * t103, 0, t119 + 0.2e1 * (-t65 * t73 * t79 + (-t65 * t130 - t79 * t133) * t81) * t81, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:55:21
	% EndTime: 2019-12-05 15:55:22
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (1050->45), mult. (1346->105), div. (247->12), fcn. (1511->9), ass. (0->61)
	t125 = sin(pkin(8));
	t127 = sin(qJ(2));
	t121 = t127 ^ 2;
	t128 = cos(qJ(2));
	t149 = t121 / t128 ^ 2;
	t138 = 0.1e1 + t149;
	t163 = t125 * t138;
	t162 = qJD(2) * t127;
	t116 = pkin(9) + qJ(4) + qJ(5);
	t114 = sin(t116);
	t115 = cos(t116);
	t126 = cos(pkin(8));
	t146 = t126 * t128;
	t103 = t125 * t114 + t115 * t146;
	t147 = t125 * t127;
	t109 = atan2(-t147, -t128);
	t106 = sin(t109);
	t107 = cos(t109);
	t95 = -t106 * t147 - t107 * t128;
	t92 = 0.1e1 / t95;
	t117 = t125 ^ 2;
	t112 = t117 * t149 + 0.1e1;
	t110 = 0.1e1 / t112;
	t97 = t110 * t163;
	t140 = t125 * t97 - 0.1e1;
	t150 = t107 * t127;
	t151 = t106 * t128;
	t153 = t125 - t97;
	t82 = -t140 * t150 - t153 * t151;
	t160 = 0.2e1 * t82;
	t99 = 0.1e1 / t103;
	t100 = 0.1e1 / t103 ^ 2;
	t93 = 0.1e1 / t95 ^ 2;
	t118 = t126 ^ 2;
	t145 = qJD(2) * t128;
	t154 = t127 * t93;
	t143 = t107 * t147;
	t96 = qJD(2) * t97;
	t81 = (-t143 + t151) * t96 + (-t125 * t151 + t150) * qJD(2);
	t157 = t81 * t92 * t93;
	t89 = t118 * t121 * t93 + 0.1e1;
	t159 = (-t121 * t157 + t145 * t154) * t118 / t89 ^ 2;
	t142 = t114 * t146;
	t102 = -t125 * t115 + t142;
	t152 = t100 * t102;
	t119 = qJD(4) + qJD(5);
	t139 = t126 * t162;
	t91 = -t119 * t142 + (t119 * t125 - t139) * t115;
	t155 = t99 * t100 * t91;
	t98 = t102 ^ 2;
	t86 = t98 * t100 + 0.1e1;
	t90 = -t103 * t119 + t114 * t139;
	t158 = (-t90 * t152 - t98 * t155) / t86 ^ 2;
	t87 = 0.1e1 / t89;
	t156 = t87 * t93;
	t144 = -0.2e1 * t158;
	t137 = -t114 * t99 + t115 * t152;
	t84 = 0.1e1 / t86;
	t83 = 0.2e1 * (t110 - t138 / t112 ^ 2 * t117) / t128 * t162 * t163;
	t78 = t144 + 0.2e1 * (-t100 * t84 * t90 + (-t100 * t158 - t84 * t155) * t102) * t102;
	t1 = [0, t83, 0, 0, 0; 0, ((-0.2e1 * t92 * t159 + (-qJD(2) * t82 - t81) * t156) * t128 + (t93 * t159 * t160 + (t83 * t93 * t143 + t157 * t160 - qJD(2) * t92 - (t153 * qJD(2) + t140 * t96) * t106 * t154) * t87 - (t106 * t83 + (qJD(2) + (-0.2e1 * t125 + t97) * t96) * t107) * t128 * t156) * t127) * t126, 0, 0, 0; 0, (t137 * t84 * t145 + (t137 * t144 + ((-0.2e1 * t102 * t155 - t119 * t99) * t115 + (-t115 * t90 + (-t102 * t119 + t91) * t114) * t100) * t84) * t127) * t126, 0, t78, t78;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end