% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S4PRRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4PRRP5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRP5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_jacobiaD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:19:37
	% EndTime: 2019-12-29 12:19:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:19:38
	% EndTime: 2019-12-29 12:19:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:19:37
	% EndTime: 2019-12-29 12:19:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:19:38
	% EndTime: 2019-12-29 12:19:38
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (429->41), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->59)
	t97 = sin(qJ(2));
	t90 = t97 ^ 2;
	t99 = cos(qJ(2));
	t123 = t90 / t99 ^ 2;
	t111 = 0.1e1 + t123;
	t94 = sin(pkin(6));
	t135 = t111 * t94;
	t134 = qJD(2) * t97;
	t95 = cos(pkin(6));
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
	t1 = [0, t56, 0, 0; 0, ((-0.2e1 * t60 * t131 + (-qJD(2) * t54 - t53) * t128) * t99 + (t61 * t131 * t132 + (t56 * t61 * t115 + t130 * t132 - qJD(2) * t60 - (-t119 * qJD(2) + t112 * t64) * t78 * t127) * t57 - (t56 * t78 + (-t112 * qJD(2) + t119 * t64) * t79) * t99 * t128) * t97) * t95, 0, 0; 0, (t110 * t97 * t116 + (t110 * t118 + ((-qJD(3) * t73 - 0.2e1 * t76 * t126) * t98 + (-t69 * t98 + (t70 - t117) * t96) * t74) * t97) * t65) * t95, t116 + 0.2e1 * (-t65 * t69 * t74 + (-t65 * t126 - t74 * t129) * t76) * t76, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:19:38
	% EndTime: 2019-12-29 12:19:38
	% DurationCPUTime: 0.66s
	% Computational Cost: add. (429->41), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->59)
	t104 = cos(qJ(2));
	t102 = sin(qJ(2));
	t95 = t102 ^ 2;
	t129 = t95 / t104 ^ 2;
	t116 = 0.1e1 + t129;
	t99 = sin(pkin(6));
	t140 = t116 * t99;
	t139 = qJD(2) * t102;
	t101 = sin(qJ(3));
	t103 = cos(qJ(3));
	t100 = cos(pkin(6));
	t121 = t100 * t104;
	t82 = t99 * t101 + t103 * t121;
	t123 = t99 * t102;
	t85 = atan2(-t123, -t104);
	t83 = sin(t85);
	t84 = cos(t85);
	t68 = -t84 * t104 - t83 * t123;
	t65 = 0.1e1 / t68;
	t78 = 0.1e1 / t82;
	t92 = t99 ^ 2;
	t88 = t92 * t129 + 0.1e1;
	t86 = 0.1e1 / t88;
	t72 = t86 * t140;
	t117 = t72 * t99 - 0.1e1;
	t125 = t84 * t102;
	t126 = t104 * t83;
	t128 = t72 - t99;
	t59 = -t117 * t125 + t128 * t126;
	t137 = 0.2e1 * t59;
	t66 = 0.1e1 / t68 ^ 2;
	t79 = 0.1e1 / t82 ^ 2;
	t120 = qJD(2) * t104;
	t127 = t102 * t66;
	t118 = t84 * t123;
	t69 = qJD(2) * t72;
	t58 = (-t118 + t126) * t69 + (-t99 * t126 + t125) * qJD(2);
	t135 = t58 * t65 * t66;
	t93 = t100 ^ 2;
	t64 = t93 * t95 * t66 + 0.1e1;
	t136 = (t120 * t127 - t95 * t135) * t93 / t64 ^ 2;
	t81 = t101 * t121 - t99 * t103;
	t130 = t79 * t81;
	t114 = t100 * t139;
	t122 = qJD(3) * t81;
	t75 = -t103 * t114 - t122;
	t131 = t75 * t78 * t79;
	t77 = t81 ^ 2;
	t73 = t77 * t79 + 0.1e1;
	t74 = -t82 * qJD(3) + t101 * t114;
	t134 = (-t74 * t130 - t77 * t131) / t73 ^ 2;
	t62 = 0.1e1 / t64;
	t133 = t62 * t66;
	t132 = t74 * t79;
	t119 = -0.2e1 * t134;
	t113 = -t101 * t78 + t103 * t130;
	t70 = 0.1e1 / t73;
	t61 = 0.2e1 * (t86 - t116 / t88 ^ 2 * t92) / t104 * t139 * t140;
	t1 = [0, t61, 0, 0; 0, ((-0.2e1 * t65 * t136 + (-qJD(2) * t59 - t58) * t133) * t104 + (t66 * t136 * t137 + (t61 * t66 * t118 + t135 * t137 - qJD(2) * t65 - (-t128 * qJD(2) + t117 * t69) * t83 * t127) * t62 - (t61 * t83 + (-t117 * qJD(2) + t128 * t69) * t84) * t104 * t133) * t102) * t100, 0, 0; 0, (t113 * t70 * t120 + (t113 * t119 + ((t75 - t122) * t79 * t101 + (-qJD(3) * t78 - 0.2e1 * t81 * t131 - t132) * t103) * t70) * t102) * t100, t119 + 0.2e1 * (-t70 * t132 + (-t70 * t131 - t79 * t134) * t81) * t81, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end