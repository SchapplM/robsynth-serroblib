% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPPR2
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
%   Wie in S5PRPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:22
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRPPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:30
	% EndTime: 2019-10-24 10:22:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:30
	% EndTime: 2019-10-24 10:22:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:31
	% EndTime: 2019-10-24 10:22:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:31
	% EndTime: 2019-10-24 10:22:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:31
	% EndTime: 2019-10-24 10:22:31
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (825->35), mult. (901->88), div. (202->12), fcn. (1017->9), ass. (0->53)
	t85 = qJ(2) + pkin(8);
	t81 = sin(t85);
	t77 = t81 ^ 2;
	t82 = cos(t85);
	t108 = t77 / t82 ^ 2;
	t100 = 0.1e1 + t108;
	t87 = sin(pkin(7));
	t118 = t100 * t87;
	t89 = cos(pkin(7));
	t84 = t89 ^ 2;
	t107 = t77 * t84;
	t117 = 0.2e1 * t107;
	t116 = qJD(2) * t81;
	t105 = t87 * t81;
	t72 = atan2(-t105, -t82);
	t70 = sin(t72);
	t71 = cos(t72);
	t57 = -t70 * t105 - t71 * t82;
	t54 = 0.1e1 / t57;
	t106 = t82 * t89;
	t86 = sin(pkin(9));
	t88 = cos(pkin(9));
	t69 = t88 * t106 + t87 * t86;
	t65 = 0.1e1 / t69;
	t55 = 0.1e1 / t57 ^ 2;
	t66 = 0.1e1 / t69 ^ 2;
	t83 = t87 ^ 2;
	t75 = t83 * t108 + 0.1e1;
	t73 = 0.1e1 / t75;
	t59 = t73 * t118;
	t109 = t70 * t82;
	t98 = -t87 * t109 + t71 * t81;
	t103 = t71 * t105;
	t99 = -t103 + t109;
	t48 = t99 * t59 + t98;
	t114 = 0.2e1 * t48;
	t111 = t55 * t82;
	t58 = qJD(2) * t59;
	t47 = t98 * qJD(2) + t99 * t58;
	t112 = t47 * t54 * t55;
	t53 = t55 * t107 + 0.1e1;
	t113 = (t111 * t116 - t77 * t112) * t84 / t53 ^ 2;
	t102 = t86 * t106;
	t68 = -t87 * t88 + t102;
	t110 = t68 * t88;
	t104 = t59 - t87;
	t101 = t59 * t87 - 0.1e1;
	t67 = t65 * t66;
	t64 = t68 ^ 2;
	t62 = t64 * t66 + 0.1e1;
	t51 = 0.1e1 / t53;
	t49 = 0.2e1 * (t73 - t100 / t75 ^ 2 * t83) / t82 * t116 * t118;
	t1 = [0, t49, 0, 0, 0; 0, ((-0.2e1 * t54 * t113 + (-qJD(2) * t48 - t47) * t55 * t51) * t82 + (t55 * t113 * t114 + (t112 * t114 - qJD(2) * t54 - (t49 * t70 + (-t101 * qJD(2) + t104 * t58) * t71) * t111 + (t49 * t103 - (-t104 * qJD(2) + t101 * t58) * t81 * t70) * t55) * t51) * t81) * t89, 0, 0, 0; 0, ((-t66 * t110 + t65 * t86) / t62 ^ 2 * (t64 * t67 * t88 - t66 * t68 * t86) * t117 + (-t65 * t102 + (t67 * t110 * t117 + (t68 * t106 - 0.2e1 * t86 * t107) * t66) * t88) / t62) * qJD(2), 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:31
	% EndTime: 2019-10-24 10:22:31
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (1158->44), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->61)
	t113 = sin(pkin(7));
	t112 = qJ(2) + pkin(8);
	t106 = sin(t112);
	t101 = t106 ^ 2;
	t108 = cos(t112);
	t134 = t101 / t108 ^ 2;
	t124 = 0.1e1 + t134;
	t150 = t113 * t124;
	t149 = qJD(2) * t106;
	t111 = pkin(9) + qJ(5);
	t105 = sin(t111);
	t107 = cos(t111);
	t114 = cos(pkin(7));
	t133 = t108 * t114;
	t90 = t113 * t105 + t107 * t133;
	t131 = t113 * t106;
	t93 = atan2(-t131, -t108);
	t91 = sin(t93);
	t92 = cos(t93);
	t76 = -t92 * t108 - t91 * t131;
	t73 = 0.1e1 / t76;
	t86 = 0.1e1 / t90;
	t109 = t113 ^ 2;
	t96 = t109 * t134 + 0.1e1;
	t94 = 0.1e1 / t96;
	t80 = t94 * t150;
	t126 = t113 * t80 - 0.1e1;
	t136 = t113 - t80;
	t137 = t92 * t106;
	t138 = t108 * t91;
	t67 = -t126 * t137 - t136 * t138;
	t147 = 0.2e1 * t67;
	t74 = 0.1e1 / t76 ^ 2;
	t87 = 0.1e1 / t90 ^ 2;
	t110 = t114 ^ 2;
	t130 = qJD(2) * t108;
	t139 = t106 * t74;
	t128 = t92 * t131;
	t77 = qJD(2) * t80;
	t66 = (-t128 + t138) * t77 + (-t113 * t138 + t137) * qJD(2);
	t145 = t66 * t73 * t74;
	t72 = t110 * t101 * t74 + 0.1e1;
	t146 = (-t101 * t145 + t130 * t139) * t110 / t72 ^ 2;
	t89 = t105 * t133 - t113 * t107;
	t140 = t87 * t89;
	t125 = t114 * t149;
	t135 = qJD(5) * t89;
	t83 = -t107 * t125 - t135;
	t141 = t83 * t86 * t87;
	t85 = t89 ^ 2;
	t81 = t85 * t87 + 0.1e1;
	t82 = -t90 * qJD(5) + t105 * t125;
	t144 = (-t82 * t140 - t85 * t141) / t81 ^ 2;
	t70 = 0.1e1 / t72;
	t143 = t70 * t74;
	t142 = t82 * t87;
	t129 = -0.2e1 * t144;
	t123 = -t105 * t86 + t107 * t140;
	t78 = 0.1e1 / t81;
	t69 = 0.2e1 * (t94 - t124 / t96 ^ 2 * t109) / t108 * t149 * t150;
	t1 = [0, t69, 0, 0, 0; 0, ((-0.2e1 * t73 * t146 + (-qJD(2) * t67 - t66) * t143) * t108 + (t74 * t146 * t147 + (t69 * t74 * t128 + t145 * t147 - qJD(2) * t73 - (t136 * qJD(2) + t126 * t77) * t91 * t139) * t70 - (t69 * t91 + (qJD(2) + (-0.2e1 * t113 + t80) * t77) * t92) * t108 * t143) * t106) * t114, 0, 0, 0; 0, (t123 * t78 * t130 + (t123 * t129 + ((t83 - t135) * t87 * t105 + (-qJD(5) * t86 - 0.2e1 * t89 * t141 - t142) * t107) * t78) * t106) * t114, 0, 0, t129 + 0.2e1 * (-t78 * t142 + (-t78 * t141 - t87 * t144) * t89) * t89;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end