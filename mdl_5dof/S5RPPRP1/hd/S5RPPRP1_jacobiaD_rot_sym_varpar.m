% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPRP1
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
%   Wie in S5RPPRP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPPRP1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
JaD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:05
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (701->40), mult. (853->109), div. (126->12), fcn. (1047->9), ass. (0->51)
	t107 = cos(pkin(8));
	t103 = 0.1e1 / t107 ^ 2;
	t106 = sin(pkin(8));
	t101 = t106 ^ 2;
	t105 = qJ(1) + pkin(7);
	t98 = sin(t105);
	t96 = t98 ^ 2;
	t92 = t96 * t101 * t103 + 0.1e1;
	t99 = cos(t105);
	t97 = t99 ^ 2;
	t127 = 0.1e1 / t92 ^ 2 * t97;
	t135 = t103 * t127;
	t90 = 0.1e1 / t92;
	t134 = (t90 - 0.1e1) * t106;
	t108 = sin(qJ(4));
	t109 = cos(qJ(4));
	t121 = t107 * t109;
	t86 = t98 * t108 + t99 * t121;
	t126 = t98 * t106;
	t89 = atan2(-t126, -t107);
	t87 = sin(t89);
	t88 = cos(t89);
	t75 = -t88 * t107 - t126 * t87;
	t70 = 0.1e1 / t75;
	t80 = 0.1e1 / t86;
	t102 = 0.1e1 / t107;
	t71 = 0.1e1 / t75 ^ 2;
	t81 = 0.1e1 / t86 ^ 2;
	t84 = t99 * t108 - t121 * t98;
	t122 = t107 * t108;
	t85 = -t98 * t109 + t122 * t99;
	t74 = qJD(1) * t84 - qJD(4) * t85;
	t129 = t74 * t80 * t81;
	t118 = t99 * t109 + t122 * t98;
	t73 = qJD(1) * t118 - t86 * qJD(4);
	t130 = t73 * t81;
	t79 = t85 ^ 2;
	t78 = t79 * t81 + 0.1e1;
	t132 = (-t79 * t129 - t85 * t130) / t78 ^ 2;
	t131 = t71 * t99;
	t128 = t84 * t85;
	t124 = qJD(1) * t98;
	t123 = t101 * t102;
	t120 = t102 * t135;
	t66 = (-t88 * t90 * t98 * t123 + t87 * t134) * t99;
	t100 = t106 * t101;
	t76 = 0.1e1 / t78;
	t72 = t70 * t71;
	t69 = t97 * t101 * t71 + 0.1e1;
	t65 = qJD(1) * t66;
	t1 = [(-t102 * t106 * t90 - 0.2e1 * t100 * t120) * t124, 0, 0, 0, 0; (0.2e1 * (t66 * t131 + t70 * t98) / t69 ^ 2 * (-t65 * t72 * t97 - t124 * t131) * t101 + ((0.2e1 * t66 * t72 * t99 + t71 * t98) * t65 + (-t99 * t70 + ((t66 + (t100 * t135 + t134) * t99 * t87) * t98 - (0.2e1 * t101 ^ 2 * t96 * t120 + (t127 + (t96 - 0.2e1 * t97) * t90) * t123) * t99 * t88) * t71) * qJD(1)) / t69) * t106, 0, 0, 0, 0; 0.2e1 * (t118 * t80 + t81 * t128) * t132 + ((-qJD(1) * t85 + qJD(4) * t84) * t80 + 0.2e1 * t128 * t129 + (t118 * t74 - (-qJD(1) * t86 + qJD(4) * t118) * t85 + t84 * t73) * t81) * t76, 0, 0, -0.2e1 * t132 + 0.2e1 * (-t76 * t130 + (-t76 * t129 - t81 * t132) * t85) * t85, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:05
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (701->40), mult. (853->108), div. (126->12), fcn. (1047->9), ass. (0->52)
	t109 = cos(pkin(8));
	t105 = 0.1e1 / t109 ^ 2;
	t108 = sin(pkin(8));
	t103 = t108 ^ 2;
	t107 = qJ(1) + pkin(7);
	t100 = sin(t107);
	t98 = t100 ^ 2;
	t94 = t98 * t103 * t105 + 0.1e1;
	t101 = cos(t107);
	t99 = t101 ^ 2;
	t131 = 0.1e1 / t94 ^ 2 * t99;
	t138 = t105 * t131;
	t92 = 0.1e1 / t94;
	t137 = (t92 - 0.1e1) * t108;
	t110 = sin(qJ(4));
	t111 = cos(qJ(4));
	t124 = t109 * t111;
	t88 = t100 * t110 + t101 * t124;
	t128 = t100 * t108;
	t91 = atan2(-t128, -t109);
	t89 = sin(t91);
	t90 = cos(t91);
	t77 = -t90 * t109 - t89 * t128;
	t72 = 0.1e1 / t77;
	t82 = 0.1e1 / t88;
	t104 = 0.1e1 / t109;
	t73 = 0.1e1 / t77 ^ 2;
	t83 = 0.1e1 / t88 ^ 2;
	t86 = -t100 * t124 + t101 * t110;
	t125 = t109 * t110;
	t87 = -t100 * t111 + t101 * t125;
	t76 = t86 * qJD(1) - t87 * qJD(4);
	t133 = t76 * t82 * t83;
	t120 = t100 * t125 + t101 * t111;
	t75 = t120 * qJD(1) - t88 * qJD(4);
	t134 = t75 * t83;
	t81 = t87 ^ 2;
	t80 = t81 * t83 + 0.1e1;
	t135 = (-t81 * t133 - t87 * t134) / t80 ^ 2;
	t132 = t86 * t87;
	t130 = t100 * t73;
	t129 = t101 * t73;
	t126 = t103 * t104;
	t123 = qJD(1) * t100;
	t122 = t104 * t138;
	t68 = (-t100 * t90 * t92 * t126 + t89 * t137) * t101;
	t102 = t108 * t103;
	t78 = 0.1e1 / t80;
	t74 = t72 * t73;
	t71 = t99 * t103 * t73 + 0.1e1;
	t67 = qJD(1) * t68;
	t1 = [(-t104 * t108 * t92 - 0.2e1 * t102 * t122) * t123, 0, 0, 0, 0; (0.2e1 * (t100 * t72 + t68 * t129) / t71 ^ 2 * (-t67 * t74 * t99 - t123 * t129) * t103 + ((0.2e1 * t101 * t68 * t74 + t130) * t67 + (t68 * t130 + (-t72 + (t102 * t138 + t137) * t89 * t130 - (0.2e1 * t103 ^ 2 * t98 * t122 + (t131 + (t98 - 0.2e1 * t99) * t92) * t126) * t73 * t90) * t101) * qJD(1)) / t71) * t108, 0, 0, 0, 0; 0.2e1 * (t120 * t82 + t132 * t83) * t135 + ((-qJD(1) * t87 + qJD(4) * t86) * t82 + 0.2e1 * t132 * t133 + (t120 * t76 - (-qJD(1) * t88 + qJD(4) * t120) * t87 + t86 * t75) * t83) * t78, 0, 0, -0.2e1 * t135 + 0.2e1 * (-t78 * t134 + (-t133 * t78 - t83 * t135) * t87) * t87, 0;];
	JaD_rot = t1;
end