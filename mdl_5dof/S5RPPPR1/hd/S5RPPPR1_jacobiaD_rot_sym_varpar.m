% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPPR1
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
%   Wie in S5RPPPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPPPR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
JaD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:23
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:23
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:22
	% EndTime: 2022-01-20 09:13:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:23
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:23
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (521->33), mult. (614->98), div. (108->12), fcn. (792->9), ass. (0->50)
	t85 = qJ(1) + pkin(7);
	t79 = cos(t85);
	t77 = t79 ^ 2;
	t87 = sin(pkin(8));
	t81 = t87 ^ 2;
	t110 = t77 * t81;
	t78 = sin(t85);
	t76 = t78 ^ 2;
	t89 = cos(pkin(8));
	t83 = 0.1e1 / t89 ^ 2;
	t75 = t76 * t81 * t83 + 0.1e1;
	t73 = 0.1e1 / t75;
	t109 = (t73 - 0.1e1) * t87;
	t104 = t78 * t87;
	t72 = atan2(-t104, -t89);
	t70 = sin(t72);
	t71 = cos(t72);
	t56 = -t70 * t104 - t71 * t89;
	t53 = 0.1e1 / t56;
	t88 = cos(pkin(9));
	t101 = t88 * t89;
	t86 = sin(pkin(9));
	t69 = t79 * t101 + t78 * t86;
	t63 = 0.1e1 / t69;
	t82 = 0.1e1 / t89;
	t54 = 0.1e1 / t56 ^ 2;
	t64 = 0.1e1 / t69 ^ 2;
	t107 = t54 * t79;
	t102 = t86 * t89;
	t68 = t79 * t102 - t78 * t88;
	t106 = t64 * t68;
	t67 = -t78 * t101 + t79 * t86;
	t105 = t67 * t68;
	t103 = t81 * t82;
	t100 = qJD(1) * t78;
	t99 = t73 * t103;
	t74 = 0.1e1 / t75 ^ 2;
	t98 = t74 * t87 * t110;
	t66 = -t78 * t102 - t79 * t88;
	t49 = (-t71 * t78 * t99 + t70 * t109) * t79;
	t84 = t82 * t83;
	t65 = t63 * t64;
	t62 = t68 ^ 2;
	t61 = t67 * qJD(1);
	t60 = t66 * qJD(1);
	t59 = t62 * t64 + 0.1e1;
	t55 = t53 * t54;
	t52 = t54 * t110 + 0.1e1;
	t48 = qJD(1) * t49;
	t1 = [(-t73 * t82 * t87 - 0.2e1 * t84 * t98) * t100, 0, 0, 0, 0; (0.2e1 * (t49 * t107 + t53 * t78) / t52 ^ 2 * (-t48 * t55 * t77 - t100 * t107) * t81 + ((0.2e1 * t49 * t55 * t79 + t54 * t78) * t48 + (-t79 * t53 + ((t49 + (t83 * t98 + t109) * t79 * t70) * t78 - (t76 * t99 + (-0.2e1 * t99 + (0.2e1 * t76 * t81 ^ 2 * t84 + t103) * t74) * t77) * t79 * t71) * t54) * qJD(1)) / t52) * t87, 0, 0, 0, 0; 0.2e1 * (t64 * t105 - t63 * t66) / t59 ^ 2 * (-t62 * t65 * t61 + t60 * t106) + (-t67 * t60 * t64 + (0.2e1 * t65 * t105 - t66 * t64) * t61 + (t69 * t106 - t68 * t63) * qJD(1)) / t59, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:23
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (895->41), mult. (853->110), div. (126->12), fcn. (1047->9), ass. (0->56)
	t115 = cos(pkin(8));
	t110 = 0.1e1 / t115 ^ 2;
	t113 = qJ(1) + pkin(7);
	t104 = sin(t113);
	t101 = t104 ^ 2;
	t114 = sin(pkin(8));
	t108 = t114 ^ 2;
	t98 = t101 * t108 * t110 + 0.1e1;
	t97 = 0.1e1 / t98 ^ 2;
	t143 = t110 * t97;
	t112 = pkin(9) + qJ(5);
	t103 = sin(t112);
	t105 = cos(t112);
	t106 = cos(t113);
	t130 = t106 * t115;
	t91 = t104 * t103 + t105 * t130;
	t132 = t104 * t114;
	t95 = atan2(-t132, -t115);
	t92 = sin(t95);
	t127 = t92 * t132;
	t93 = cos(t95);
	t83 = -t115 * t93 - t127;
	t80 = 0.1e1 / t83;
	t85 = 0.1e1 / t91;
	t109 = 0.1e1 / t115;
	t81 = 0.1e1 / t83 ^ 2;
	t86 = 0.1e1 / t91 ^ 2;
	t96 = 0.1e1 / t98;
	t142 = t96 - 0.1e1;
	t131 = t104 * t115;
	t89 = t103 * t106 - t105 * t131;
	t90 = t103 * t130 - t104 * t105;
	t76 = qJD(1) * t89 - qJD(5) * t90;
	t139 = t76 * t85 * t86;
	t124 = t103 * t131 + t105 * t106;
	t75 = t124 * qJD(1) - qJD(5) * t91;
	t140 = t75 * t86;
	t84 = t90 ^ 2;
	t79 = t84 * t86 + 0.1e1;
	t141 = (-t139 * t84 - t140 * t90) / t79 ^ 2;
	t138 = t89 * t90;
	t102 = t106 ^ 2;
	t137 = t102 * t81;
	t136 = t104 * t81;
	t135 = t106 * t81;
	t134 = t109 * t143;
	t129 = t108 * t109;
	t128 = qJD(1) * t104;
	t125 = t93 * t96 * t129;
	t71 = (t114 * t142 * t92 - t104 * t125) * t106;
	t107 = t114 * t108;
	t82 = t80 * t81;
	t77 = 0.1e1 / t79;
	t74 = t108 * t137 + 0.1e1;
	t70 = qJD(1) * t71;
	t1 = [(-0.2e1 * t102 * t107 * t134 - t109 * t114 * t96) * t128, 0, 0, 0, 0; (0.2e1 * (t104 * t80 + t135 * t71) / t74 ^ 2 * (-t102 * t70 * t82 - t128 * t135) * t108 + ((0.2e1 * t106 * t71 * t82 + t136) * t70 + (t71 * t136 + (-t80 - (-t104 * t107 * t92 * t143 + (0.2e1 * t101 * t108 ^ 2 * t134 + (-0.2e1 * t96 + t97) * t129) * t93) * t137 + (-t101 * t125 + t127 * t142) * t81) * t106) * qJD(1)) / t74) * t114, 0, 0, 0, 0; 0.2e1 * (t124 * t85 + t138 * t86) * t141 + ((-qJD(1) * t90 + qJD(5) * t89) * t85 + 0.2e1 * t138 * t139 + (t124 * t76 - (-qJD(1) * t91 + qJD(5) * t124) * t90 + t89 * t75) * t86) * t77, 0, 0, 0, -0.2e1 * t141 + 0.2e1 * (-t77 * t140 + (-t139 * t77 - t141 * t86) * t90) * t90;];
	JaD_rot = t1;
end