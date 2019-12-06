% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPPR3
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
%   Wie in S5PRPPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRPPR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:27:13
	% EndTime: 2019-12-05 15:27:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:27:13
	% EndTime: 2019-12-05 15:27:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:27:13
	% EndTime: 2019-12-05 15:27:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:27:13
	% EndTime: 2019-12-05 15:27:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:27:13
	% EndTime: 2019-12-05 15:27:14
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (752->28), mult. (758->66), div. (186->13), fcn. (854->7), ass. (0->42)
	t66 = qJ(2) + pkin(8);
	t62 = sin(t66);
	t58 = t62 ^ 2;
	t63 = cos(t66);
	t72 = t63 ^ 2;
	t87 = t58 / t72;
	t79 = 0.1e1 + t87;
	t67 = sin(pkin(7));
	t85 = t67 * t62;
	t52 = atan2(-t85, -t63);
	t50 = sin(t52);
	t51 = cos(t52);
	t46 = -t50 * t85 - t51 * t63;
	t43 = 0.1e1 / t46;
	t44 = 0.1e1 / t46 ^ 2;
	t64 = t67 ^ 2;
	t55 = t64 * t87 + 0.1e1;
	t53 = 0.1e1 / t55;
	t48 = t79 * t67 * t53;
	t88 = t50 * t63;
	t77 = t51 * t62 - t67 * t88;
	t81 = t51 * t85;
	t78 = -t81 + t88;
	t38 = t78 * t48 + t77;
	t93 = 0.2e1 * t38;
	t68 = cos(pkin(7));
	t65 = t68 ^ 2;
	t86 = t58 * t65;
	t42 = t44 * t86 + 0.1e1;
	t82 = qJD(2) * t63;
	t89 = t44 * t62;
	t47 = qJD(2) * t48;
	t37 = t77 * qJD(2) + t78 * t47;
	t91 = t37 * t43 * t44;
	t92 = (-t58 * t91 + t82 * t89) * t65 / t42 ^ 2;
	t40 = 0.1e1 / t42;
	t90 = t40 * t44;
	t83 = t48 - t67;
	t80 = t48 * t67 - 0.1e1;
	t56 = 0.1e1 + t65 * t72 / t67 ^ 2;
	t39 = 0.2e1 * (t53 - t79 / t55 ^ 2 * t64) * qJD(2) * t79 / t63 * t85;
	t1 = [0, t39, 0, 0, 0; 0, ((-0.2e1 * t43 * t92 + (-qJD(2) * t38 - t37) * t90) * t63 + (t44 * t92 * t93 + (t39 * t44 * t81 + t91 * t93 - qJD(2) * t43 - (-t83 * qJD(2) + t80 * t47) * t50 * t89) * t40 - (t39 * t50 + (-t80 * qJD(2) + t83 * t47) * t51) * t63 * t90) * t62) * t68, 0, 0, 0; 0, (-0.1e1 / t56 - 0.2e1 / t64 / t56 ^ 2 * t86) * t68 / t67 * t82, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:27:14
	% EndTime: 2019-12-05 15:27:14
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (899->44), mult. (1166->103), div. (229->12), fcn. (1315->9), ass. (0->54)
	t103 = sin(pkin(7));
	t102 = qJ(2) + pkin(8);
	t99 = cos(t102);
	t97 = t99 ^ 2;
	t98 = sin(t102);
	t124 = 0.1e1 / t98 ^ 2 * t97;
	t118 = 0.1e1 + t124;
	t135 = t103 * t118;
	t121 = qJD(2) * t99;
	t105 = sin(qJ(5));
	t106 = cos(qJ(5));
	t104 = cos(pkin(7));
	t122 = t104 * t98;
	t115 = -t103 * t105 + t106 * t122;
	t134 = t115 * qJD(5);
	t123 = t103 * t99;
	t88 = atan2(-t123, t98);
	t86 = sin(t88);
	t87 = cos(t88);
	t72 = -t86 * t123 + t87 * t98;
	t69 = 0.1e1 / t72;
	t85 = t103 * t106 + t105 * t122;
	t81 = 0.1e1 / t85;
	t70 = 0.1e1 / t72 ^ 2;
	t82 = 0.1e1 / t85 ^ 2;
	t119 = t87 * t123;
	t125 = t87 * t99;
	t126 = t86 * t98;
	t100 = t103 ^ 2;
	t91 = t100 * t124 + 0.1e1;
	t89 = 0.1e1 / t91;
	t73 = t89 * t135;
	t68 = qJD(2) * t73;
	t61 = (-t119 - t126) * t68 + (t103 * t126 + t125) * qJD(2);
	t132 = t61 * t69 * t70;
	t127 = t82 * t115;
	t117 = t104 * t121;
	t77 = t105 * t117 + t134;
	t129 = t77 * t81 * t82;
	t80 = t115 ^ 2;
	t76 = t80 * t82 + 0.1e1;
	t78 = t85 * qJD(5) - t106 * t117;
	t131 = (-t78 * t127 - t80 * t129) / t76 ^ 2;
	t130 = t70 * t98;
	t128 = t78 * t82;
	t101 = t104 ^ 2;
	t67 = t101 * t97 * t70 + 0.1e1;
	t120 = 0.2e1 * (-t121 * t130 - t97 * t132) * t101 / t67 ^ 2;
	t116 = -t105 * t127 + t106 * t81;
	t74 = 0.1e1 / t76;
	t65 = 0.1e1 / t67;
	t63 = -0.2e1 * (t89 - t118 / t91 ^ 2 * t100) / t98 * t121 * t135;
	t62 = -t73 * t126 + t125 + (-t73 * t125 + t126) * t103;
	t1 = [0, t63, 0, 0, 0; 0, ((t69 * t120 + (qJD(2) * t62 + t61) * t70 * t65) * t98 + (t62 * t70 * t120 + (0.2e1 * t62 * t132 - qJD(2) * t69 - (-t63 * t86 + (-qJD(2) + (0.2e1 * t103 - t73) * t68) * t87) * t130 + (t63 * t119 - ((t103 * t73 - 0.1e1) * t68 + (t103 - t73) * qJD(2)) * t99 * t86) * t70) * t65) * t99) * t104, 0, 0, 0; 0, (0.2e1 * t116 * t99 * t131 + (t116 * t98 * qJD(2) + ((t77 + t134) * t82 * t106 + (qJD(5) * t81 - 0.2e1 * t115 * t129 - t128) * t105) * t99) * t74) * t104, 0, 0, -0.2e1 * t131 - 0.2e1 * (t74 * t128 - (-t74 * t129 - t82 * t131) * t115) * t115;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end