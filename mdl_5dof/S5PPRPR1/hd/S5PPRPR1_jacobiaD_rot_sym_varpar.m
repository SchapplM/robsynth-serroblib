% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRPR1
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
%   Wie in S5PPRPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PPRPR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:49
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (825->35), mult. (901->88), div. (202->12), fcn. (1017->9), ass. (0->53)
	t86 = sin(pkin(7));
	t84 = pkin(8) + qJ(3);
	t80 = sin(t84);
	t76 = t80 ^ 2;
	t81 = cos(t84);
	t107 = t76 / t81 ^ 2;
	t99 = 0.1e1 + t107;
	t117 = t86 * t99;
	t88 = cos(pkin(7));
	t83 = t88 ^ 2;
	t106 = t76 * t83;
	t116 = 0.2e1 * t106;
	t115 = qJD(3) * t80;
	t104 = t86 * t80;
	t71 = atan2(-t104, -t81);
	t69 = sin(t71);
	t70 = cos(t71);
	t56 = -t69 * t104 - t70 * t81;
	t53 = 0.1e1 / t56;
	t105 = t81 * t88;
	t85 = sin(pkin(9));
	t87 = cos(pkin(9));
	t68 = t87 * t105 + t86 * t85;
	t64 = 0.1e1 / t68;
	t54 = 0.1e1 / t56 ^ 2;
	t65 = 0.1e1 / t68 ^ 2;
	t82 = t86 ^ 2;
	t74 = t82 * t107 + 0.1e1;
	t72 = 0.1e1 / t74;
	t58 = t72 * t117;
	t108 = t69 * t81;
	t97 = -t86 * t108 + t70 * t80;
	t102 = t70 * t104;
	t98 = -t102 + t108;
	t47 = t98 * t58 + t97;
	t113 = 0.2e1 * t47;
	t110 = t54 * t81;
	t57 = qJD(3) * t58;
	t46 = t97 * qJD(3) + t98 * t57;
	t111 = t46 * t53 * t54;
	t52 = t54 * t106 + 0.1e1;
	t112 = (t110 * t115 - t76 * t111) * t83 / t52 ^ 2;
	t101 = t85 * t105;
	t67 = -t86 * t87 + t101;
	t109 = t67 * t87;
	t103 = t58 - t86;
	t100 = t58 * t86 - 0.1e1;
	t66 = t64 * t65;
	t63 = t67 ^ 2;
	t61 = t63 * t65 + 0.1e1;
	t50 = 0.1e1 / t52;
	t48 = 0.2e1 * (t72 - t99 / t74 ^ 2 * t82) / t81 * t115 * t117;
	t1 = [0, 0, t48, 0, 0; 0, 0, ((-0.2e1 * t53 * t112 + (-qJD(3) * t47 - t46) * t54 * t50) * t81 + (t54 * t112 * t113 + (t111 * t113 - qJD(3) * t53 - (t48 * t69 + (-t100 * qJD(3) + t103 * t57) * t70) * t110 + (t48 * t102 - (-t103 * qJD(3) + t100 * t57) * t80 * t69) * t54) * t50) * t80) * t88, 0, 0; 0, 0, ((-t65 * t109 + t64 * t85) / t61 ^ 2 * (t63 * t66 * t87 - t65 * t67 * t85) * t116 + (-t64 * t101 + (t66 * t109 * t116 + (t67 * t105 - 0.2e1 * t85 * t106) * t65) * t87) / t61) * qJD(3), 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:49
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (1158->44), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->61)
	t114 = sin(pkin(7));
	t113 = pkin(8) + qJ(3);
	t107 = sin(t113);
	t102 = t107 ^ 2;
	t109 = cos(t113);
	t135 = t102 / t109 ^ 2;
	t125 = 0.1e1 + t135;
	t151 = t114 * t125;
	t150 = qJD(3) * t107;
	t112 = pkin(9) + qJ(5);
	t106 = sin(t112);
	t108 = cos(t112);
	t115 = cos(pkin(7));
	t134 = t109 * t115;
	t91 = t114 * t106 + t108 * t134;
	t132 = t114 * t107;
	t94 = atan2(-t132, -t109);
	t92 = sin(t94);
	t93 = cos(t94);
	t77 = -t93 * t109 - t92 * t132;
	t74 = 0.1e1 / t77;
	t87 = 0.1e1 / t91;
	t110 = t114 ^ 2;
	t97 = t110 * t135 + 0.1e1;
	t95 = 0.1e1 / t97;
	t81 = t95 * t151;
	t127 = t114 * t81 - 0.1e1;
	t137 = t114 - t81;
	t138 = t93 * t107;
	t139 = t109 * t92;
	t68 = -t127 * t138 - t137 * t139;
	t148 = 0.2e1 * t68;
	t75 = 0.1e1 / t77 ^ 2;
	t88 = 0.1e1 / t91 ^ 2;
	t111 = t115 ^ 2;
	t131 = qJD(3) * t109;
	t140 = t107 * t75;
	t129 = t93 * t132;
	t78 = qJD(3) * t81;
	t67 = (-t129 + t139) * t78 + (-t114 * t139 + t138) * qJD(3);
	t146 = t67 * t74 * t75;
	t73 = t111 * t102 * t75 + 0.1e1;
	t147 = (-t102 * t146 + t131 * t140) * t111 / t73 ^ 2;
	t90 = t106 * t134 - t114 * t108;
	t141 = t88 * t90;
	t126 = t115 * t150;
	t136 = qJD(5) * t90;
	t84 = -t108 * t126 - t136;
	t142 = t84 * t87 * t88;
	t86 = t90 ^ 2;
	t82 = t86 * t88 + 0.1e1;
	t83 = -t91 * qJD(5) + t106 * t126;
	t145 = (-t83 * t141 - t86 * t142) / t82 ^ 2;
	t71 = 0.1e1 / t73;
	t144 = t71 * t75;
	t143 = t83 * t88;
	t130 = -0.2e1 * t145;
	t124 = -t106 * t87 + t108 * t141;
	t79 = 0.1e1 / t82;
	t70 = 0.2e1 * (t95 - t125 / t97 ^ 2 * t110) / t109 * t150 * t151;
	t1 = [0, 0, t70, 0, 0; 0, 0, ((-0.2e1 * t74 * t147 + (-qJD(3) * t68 - t67) * t144) * t109 + (t75 * t147 * t148 + (t70 * t75 * t129 + t146 * t148 - qJD(3) * t74 - (t137 * qJD(3) + t127 * t78) * t92 * t140) * t71 - (t70 * t92 + (qJD(3) + (-0.2e1 * t114 + t81) * t78) * t93) * t109 * t144) * t107) * t115, 0, 0; 0, 0, (t124 * t79 * t131 + (t124 * t130 + ((t84 - t136) * t88 * t106 + (-qJD(5) * t87 - 0.2e1 * t90 * t142 - t143) * t108) * t79) * t107) * t115, 0, t130 + 0.2e1 * (-t79 * t143 + (-t79 * t142 - t88 * t145) * t90) * t90;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end