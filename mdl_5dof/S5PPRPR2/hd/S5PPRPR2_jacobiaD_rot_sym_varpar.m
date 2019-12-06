% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRPR2
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
%   Wie in S5PPRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PPRPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:39
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (752->28), mult. (758->66), div. (186->13), fcn. (854->7), ass. (0->42)
	t65 = pkin(8) + qJ(3);
	t61 = sin(t65);
	t57 = t61 ^ 2;
	t62 = cos(t65);
	t71 = t62 ^ 2;
	t86 = t57 / t71;
	t78 = 0.1e1 + t86;
	t66 = sin(pkin(7));
	t84 = t66 * t61;
	t51 = atan2(-t84, -t62);
	t49 = sin(t51);
	t50 = cos(t51);
	t45 = -t49 * t84 - t50 * t62;
	t42 = 0.1e1 / t45;
	t43 = 0.1e1 / t45 ^ 2;
	t63 = t66 ^ 2;
	t54 = t63 * t86 + 0.1e1;
	t52 = 0.1e1 / t54;
	t47 = t78 * t66 * t52;
	t87 = t49 * t62;
	t76 = t50 * t61 - t66 * t87;
	t80 = t50 * t84;
	t77 = -t80 + t87;
	t37 = t77 * t47 + t76;
	t92 = 0.2e1 * t37;
	t67 = cos(pkin(7));
	t64 = t67 ^ 2;
	t85 = t57 * t64;
	t41 = t43 * t85 + 0.1e1;
	t81 = qJD(3) * t62;
	t88 = t43 * t61;
	t46 = qJD(3) * t47;
	t36 = t76 * qJD(3) + t77 * t46;
	t90 = t36 * t42 * t43;
	t91 = (-t57 * t90 + t81 * t88) * t64 / t41 ^ 2;
	t39 = 0.1e1 / t41;
	t89 = t39 * t43;
	t82 = t47 - t66;
	t79 = t47 * t66 - 0.1e1;
	t55 = 0.1e1 + t64 * t71 / t66 ^ 2;
	t38 = 0.2e1 * (t52 - t78 / t54 ^ 2 * t63) * qJD(3) * t78 / t62 * t84;
	t1 = [0, 0, t38, 0, 0; 0, 0, ((-0.2e1 * t42 * t91 + (-qJD(3) * t37 - t36) * t89) * t62 + (t43 * t91 * t92 + (t38 * t43 * t80 + t90 * t92 - qJD(3) * t42 - (-t82 * qJD(3) + t79 * t46) * t49 * t88) * t39 - (t38 * t49 + (-t79 * qJD(3) + t82 * t46) * t50) * t62 * t89) * t61) * t67, 0, 0; 0, 0, (-0.1e1 / t55 - 0.2e1 / t63 / t55 ^ 2 * t85) * t67 / t66 * t81, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:39
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (899->44), mult. (1166->103), div. (229->12), fcn. (1315->9), ass. (0->54)
	t102 = sin(pkin(7));
	t101 = pkin(8) + qJ(3);
	t98 = cos(t101);
	t96 = t98 ^ 2;
	t97 = sin(t101);
	t123 = 0.1e1 / t97 ^ 2 * t96;
	t117 = 0.1e1 + t123;
	t134 = t102 * t117;
	t120 = qJD(3) * t98;
	t104 = sin(qJ(5));
	t105 = cos(qJ(5));
	t103 = cos(pkin(7));
	t121 = t103 * t97;
	t114 = -t102 * t104 + t105 * t121;
	t133 = t114 * qJD(5);
	t122 = t102 * t98;
	t87 = atan2(-t122, t97);
	t85 = sin(t87);
	t86 = cos(t87);
	t71 = -t85 * t122 + t86 * t97;
	t68 = 0.1e1 / t71;
	t84 = t102 * t105 + t104 * t121;
	t80 = 0.1e1 / t84;
	t69 = 0.1e1 / t71 ^ 2;
	t81 = 0.1e1 / t84 ^ 2;
	t118 = t86 * t122;
	t124 = t86 * t98;
	t125 = t85 * t97;
	t99 = t102 ^ 2;
	t90 = t99 * t123 + 0.1e1;
	t88 = 0.1e1 / t90;
	t72 = t88 * t134;
	t67 = qJD(3) * t72;
	t60 = (-t118 - t125) * t67 + (t102 * t125 + t124) * qJD(3);
	t131 = t60 * t68 * t69;
	t126 = t81 * t114;
	t116 = t103 * t120;
	t76 = t104 * t116 + t133;
	t128 = t76 * t80 * t81;
	t79 = t114 ^ 2;
	t75 = t79 * t81 + 0.1e1;
	t77 = t84 * qJD(5) - t105 * t116;
	t130 = (-t77 * t126 - t79 * t128) / t75 ^ 2;
	t129 = t69 * t97;
	t127 = t77 * t81;
	t100 = t103 ^ 2;
	t66 = t100 * t96 * t69 + 0.1e1;
	t119 = 0.2e1 * (-t120 * t129 - t96 * t131) * t100 / t66 ^ 2;
	t115 = -t104 * t126 + t105 * t80;
	t73 = 0.1e1 / t75;
	t64 = 0.1e1 / t66;
	t62 = -0.2e1 * (t88 - t117 / t90 ^ 2 * t99) / t97 * t120 * t134;
	t61 = -t72 * t125 + t124 + (-t72 * t124 + t125) * t102;
	t1 = [0, 0, t62, 0, 0; 0, 0, ((t68 * t119 + (qJD(3) * t61 + t60) * t69 * t64) * t97 + (t61 * t69 * t119 + (0.2e1 * t61 * t131 - qJD(3) * t68 - (-t62 * t85 + (-qJD(3) + (0.2e1 * t102 - t72) * t67) * t86) * t129 + (t62 * t118 - ((t102 * t72 - 0.1e1) * t67 + (t102 - t72) * qJD(3)) * t98 * t85) * t69) * t64) * t98) * t103, 0, 0; 0, 0, (0.2e1 * t115 * t98 * t130 + (t115 * t97 * qJD(3) + ((t76 + t133) * t81 * t105 + (qJD(5) * t80 - 0.2e1 * t114 * t128 - t127) * t104) * t98) * t73) * t103, 0, -0.2e1 * t130 - 0.2e1 * (t73 * t127 - (-t73 * t128 - t81 * t130) * t114) * t114;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end