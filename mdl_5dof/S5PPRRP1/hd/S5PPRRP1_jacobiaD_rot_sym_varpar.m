% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRRP1
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
%   Wie in S5PPRRP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PPRRP1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:07:31
	% EndTime: 2019-12-05 15:07:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:07:31
	% EndTime: 2019-12-05 15:07:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:07:31
	% EndTime: 2019-12-05 15:07:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:07:31
	% EndTime: 2019-12-05 15:07:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:07:31
	% EndTime: 2019-12-05 15:07:32
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (973->43), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->60)
	t104 = sin(pkin(7));
	t103 = pkin(8) + qJ(3);
	t100 = cos(t103);
	t99 = sin(t103);
	t95 = t99 ^ 2;
	t130 = t95 / t100 ^ 2;
	t120 = 0.1e1 + t130;
	t143 = t104 * t120;
	t142 = qJD(3) * t99;
	t106 = sin(qJ(4));
	t107 = cos(qJ(4));
	t105 = cos(pkin(7));
	t125 = t100 * t105;
	t84 = t104 * t106 + t107 * t125;
	t128 = t104 * t99;
	t87 = atan2(-t128, -t100);
	t85 = sin(t87);
	t86 = cos(t87);
	t70 = -t86 * t100 - t85 * t128;
	t67 = 0.1e1 / t70;
	t80 = 0.1e1 / t84;
	t101 = t104 ^ 2;
	t90 = t101 * t130 + 0.1e1;
	t88 = 0.1e1 / t90;
	t72 = t88 * t143;
	t117 = t104 * t72 - 0.1e1;
	t127 = t104 - t72;
	t129 = t100 * t85;
	t131 = t86 * t99;
	t61 = -t117 * t131 - t127 * t129;
	t140 = 0.2e1 * t61;
	t68 = 0.1e1 / t70 ^ 2;
	t81 = 0.1e1 / t84 ^ 2;
	t102 = t105 ^ 2;
	t123 = qJD(3) * t100;
	t135 = t68 * t99;
	t121 = t86 * t128;
	t71 = qJD(3) * t72;
	t60 = (-t121 + t129) * t71 + (-t104 * t129 + t131) * qJD(3);
	t138 = t60 * t67 * t68;
	t66 = t102 * t95 * t68 + 0.1e1;
	t139 = (t123 * t135 - t95 * t138) * t102 / t66 ^ 2;
	t83 = -t104 * t107 + t106 * t125;
	t132 = t81 * t83;
	t119 = t105 * t142;
	t126 = qJD(4) * t83;
	t77 = -t107 * t119 - t126;
	t133 = t77 * t80 * t81;
	t79 = t83 ^ 2;
	t75 = t79 * t81 + 0.1e1;
	t76 = -t84 * qJD(4) + t106 * t119;
	t137 = (-t76 * t132 - t79 * t133) / t75 ^ 2;
	t64 = 0.1e1 / t66;
	t136 = t64 * t68;
	t134 = t76 * t81;
	t122 = -0.2e1 * t137;
	t116 = -t106 * t80 + t107 * t132;
	t73 = 0.1e1 / t75;
	t62 = 0.2e1 * (t88 - t120 / t90 ^ 2 * t101) / t100 * t142 * t143;
	t1 = [0, 0, t62, 0, 0; 0, 0, ((-0.2e1 * t67 * t139 + (-qJD(3) * t61 - t60) * t136) * t100 + (t68 * t139 * t140 + (t62 * t68 * t121 + t138 * t140 - qJD(3) * t67 - (t127 * qJD(3) + t117 * t71) * t85 * t135) * t64 - (t62 * t85 + (qJD(3) + (-0.2e1 * t104 + t72) * t71) * t86) * t100 * t136) * t99) * t105, 0, 0; 0, 0, (t116 * t99 * t122 + (t116 * t123 + ((t77 - t126) * t81 * t106 + (-qJD(4) * t80 - 0.2e1 * t83 * t133 - t134) * t107) * t99) * t73) * t105, t122 + 0.2e1 * (-t73 * t134 + (-t73 * t133 - t81 * t137) * t83) * t83, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:07:32
	% EndTime: 2019-12-05 15:07:32
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (973->43), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->60)
	t109 = sin(pkin(7));
	t108 = pkin(8) + qJ(3);
	t104 = sin(t108);
	t100 = t104 ^ 2;
	t105 = cos(t108);
	t132 = t100 / t105 ^ 2;
	t122 = 0.1e1 + t132;
	t148 = t109 * t122;
	t147 = qJD(3) * t104;
	t111 = sin(qJ(4));
	t112 = cos(qJ(4));
	t110 = cos(pkin(7));
	t131 = t105 * t110;
	t89 = t109 * t111 + t112 * t131;
	t130 = t109 * t104;
	t92 = atan2(-t130, -t105);
	t90 = sin(t92);
	t91 = cos(t92);
	t75 = -t91 * t105 - t90 * t130;
	t72 = 0.1e1 / t75;
	t85 = 0.1e1 / t89;
	t106 = t109 ^ 2;
	t95 = t106 * t132 + 0.1e1;
	t93 = 0.1e1 / t95;
	t77 = t93 * t148;
	t124 = t109 * t77 - 0.1e1;
	t134 = t109 - t77;
	t135 = t91 * t104;
	t136 = t105 * t90;
	t66 = -t124 * t135 - t134 * t136;
	t145 = 0.2e1 * t66;
	t73 = 0.1e1 / t75 ^ 2;
	t86 = 0.1e1 / t89 ^ 2;
	t107 = t110 ^ 2;
	t128 = qJD(3) * t105;
	t137 = t104 * t73;
	t126 = t91 * t130;
	t76 = qJD(3) * t77;
	t65 = (-t126 + t136) * t76 + (-t109 * t136 + t135) * qJD(3);
	t143 = t65 * t72 * t73;
	t71 = t107 * t100 * t73 + 0.1e1;
	t144 = (-t100 * t143 + t128 * t137) * t107 / t71 ^ 2;
	t88 = -t109 * t112 + t111 * t131;
	t138 = t86 * t88;
	t123 = t110 * t147;
	t133 = qJD(4) * t88;
	t82 = -t112 * t123 - t133;
	t139 = t82 * t85 * t86;
	t84 = t88 ^ 2;
	t80 = t84 * t86 + 0.1e1;
	t81 = -t89 * qJD(4) + t111 * t123;
	t142 = (-t81 * t138 - t84 * t139) / t80 ^ 2;
	t69 = 0.1e1 / t71;
	t141 = t69 * t73;
	t140 = t81 * t86;
	t127 = -0.2e1 * t142;
	t121 = -t111 * t85 + t112 * t138;
	t78 = 0.1e1 / t80;
	t67 = 0.2e1 * (t93 - t122 / t95 ^ 2 * t106) / t105 * t147 * t148;
	t1 = [0, 0, t67, 0, 0; 0, 0, ((-0.2e1 * t72 * t144 + (-qJD(3) * t66 - t65) * t141) * t105 + (t73 * t144 * t145 + (t67 * t73 * t126 + t143 * t145 - qJD(3) * t72 - (t134 * qJD(3) + t124 * t76) * t90 * t137) * t69 - (t67 * t90 + (qJD(3) + (-0.2e1 * t109 + t77) * t76) * t91) * t105 * t141) * t104) * t110, 0, 0; 0, 0, (t121 * t78 * t128 + (t121 * t127 + ((t82 - t133) * t86 * t111 + (-qJD(4) * t85 - 0.2e1 * t88 * t139 - t140) * t112) * t78) * t104) * t110, t127 + 0.2e1 * (-t78 * t140 + (-t78 * t139 - t86 * t142) * t88) * t88, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end