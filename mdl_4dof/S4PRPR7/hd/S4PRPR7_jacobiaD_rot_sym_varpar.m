% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4PRPR7
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
%   Wie in S4PRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4PRPR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRPR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_jacobiaD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:26:01
	% EndTime: 2019-12-31 16:26:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:26:01
	% EndTime: 2019-12-31 16:26:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:26:01
	% EndTime: 2019-12-31 16:26:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:26:01
	% EndTime: 2019-12-31 16:26:01
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (307->27), mult. (758->66), div. (186->13), fcn. (854->7), ass. (0->41)
	t62 = sin(qJ(2));
	t56 = t62 ^ 2;
	t63 = cos(qJ(2));
	t69 = t63 ^ 2;
	t81 = t56 / t69;
	t74 = 0.1e1 + t81;
	t60 = sin(pkin(6));
	t80 = t60 * t62;
	t48 = atan2(-t80, -t63);
	t46 = sin(t48);
	t47 = cos(t48);
	t42 = -t46 * t80 - t47 * t63;
	t39 = 0.1e1 / t42;
	t40 = 0.1e1 / t42 ^ 2;
	t53 = t60 ^ 2;
	t52 = t53 * t81 + 0.1e1;
	t49 = 0.1e1 / t52;
	t44 = t74 * t60 * t49;
	t83 = t46 * t63;
	t72 = t47 * t62 - t60 * t83;
	t76 = t47 * t80;
	t73 = -t76 + t83;
	t34 = t73 * t44 + t72;
	t88 = 0.2e1 * t34;
	t61 = cos(pkin(6));
	t54 = t61 ^ 2;
	t82 = t54 * t56;
	t38 = t40 * t82 + 0.1e1;
	t77 = qJD(2) * t63;
	t84 = t40 * t62;
	t43 = qJD(2) * t44;
	t33 = t72 * qJD(2) + t73 * t43;
	t86 = t33 * t39 * t40;
	t87 = (-t56 * t86 + t77 * t84) * t54 / t38 ^ 2;
	t36 = 0.1e1 / t38;
	t85 = t36 * t40;
	t78 = t44 - t60;
	t75 = t44 * t60 - 0.1e1;
	t51 = 0.1e1 + t54 * t69 / t60 ^ 2;
	t35 = 0.2e1 * (t49 - t74 / t52 ^ 2 * t53) * qJD(2) * t74 / t63 * t80;
	t1 = [0, t35, 0, 0; 0, ((-0.2e1 * t39 * t87 + (-qJD(2) * t34 - t33) * t85) * t63 + (t40 * t87 * t88 + (t35 * t40 * t76 + t86 * t88 - qJD(2) * t39 - (-t78 * qJD(2) + t75 * t43) * t46 * t84) * t36 - (t35 * t46 + (-t75 * qJD(2) + t78 * t43) * t47) * t63 * t85) * t62) * t61, 0, 0; 0, (-0.1e1 / t51 - 0.2e1 / t53 / t51 ^ 2 * t82) * t61 / t60 * t77, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:26:01
	% EndTime: 2019-12-31 16:26:02
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (355->40), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->55)
	t97 = cos(qJ(2));
	t91 = t97 ^ 2;
	t95 = sin(qJ(2));
	t119 = 0.1e1 / t95 ^ 2 * t91;
	t110 = 0.1e1 + t119;
	t92 = sin(pkin(6));
	t128 = t110 * t92;
	t115 = qJD(2) * t97;
	t93 = cos(pkin(6));
	t117 = t93 * t95;
	t94 = sin(qJ(4));
	t96 = cos(qJ(4));
	t106 = t96 * t117 - t92 * t94;
	t127 = t106 * qJD(4);
	t118 = t92 * t97;
	t80 = atan2(-t118, t95);
	t78 = sin(t80);
	t79 = cos(t80);
	t64 = -t78 * t118 + t79 * t95;
	t61 = 0.1e1 / t64;
	t77 = t94 * t117 + t92 * t96;
	t73 = 0.1e1 / t77;
	t62 = 0.1e1 / t64 ^ 2;
	t74 = 0.1e1 / t77 ^ 2;
	t120 = t78 * t95;
	t107 = t92 * t120 + t79 * t97;
	t113 = t79 * t118;
	t108 = -t113 - t120;
	t85 = t92 ^ 2;
	t83 = t85 * t119 + 0.1e1;
	t81 = 0.1e1 / t83;
	t67 = t81 * t128;
	t60 = qJD(2) * t67;
	t53 = t107 * qJD(2) + t108 * t60;
	t125 = t53 * t61 * t62;
	t121 = t74 * t106;
	t112 = t93 * t115;
	t69 = t94 * t112 + t127;
	t122 = t69 * t73 * t74;
	t72 = t106 ^ 2;
	t68 = t72 * t74 + 0.1e1;
	t70 = t77 * qJD(4) - t96 * t112;
	t124 = (-t70 * t121 - t72 * t122) / t68 ^ 2;
	t123 = t62 * t95;
	t116 = -t67 + t92;
	t86 = t93 ^ 2;
	t59 = t86 * t91 * t62 + 0.1e1;
	t114 = 0.2e1 * (-t115 * t123 - t91 * t125) * t86 / t59 ^ 2;
	t111 = t67 * t92 - 0.1e1;
	t109 = -t94 * t121 + t73 * t96;
	t65 = 0.1e1 / t68;
	t57 = 0.1e1 / t59;
	t56 = -0.2e1 * (t81 - t110 / t83 ^ 2 * t85) / t95 * t115 * t128;
	t54 = t108 * t67 + t107;
	t1 = [0, t56, 0, 0; 0, ((t61 * t114 + (qJD(2) * t54 + t53) * t62 * t57) * t95 + (t54 * t62 * t114 + (0.2e1 * t54 * t125 - qJD(2) * t61 - (-t56 * t78 + (t111 * qJD(2) + t116 * t60) * t79) * t123 + (t56 * t113 - (t116 * qJD(2) + t111 * t60) * t97 * t78) * t62) * t57) * t97) * t93, 0, 0; 0, (0.2e1 * t109 * t97 * t124 + (t109 * t95 * qJD(2) + ((qJD(4) * t73 - 0.2e1 * t106 * t122) * t94 + (-t70 * t94 + (t69 + t127) * t96) * t74) * t97) * t65) * t93, 0, -0.2e1 * t124 - 0.2e1 * (t65 * t70 * t74 - (-t65 * t122 - t74 * t124) * t106) * t106;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end