% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPRR7
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
%   Wie in S5PRPRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:28
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRPRR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:59
	% EndTime: 2019-10-24 10:27:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:59
	% EndTime: 2019-10-24 10:27:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:59
	% EndTime: 2019-10-24 10:27:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:59
	% EndTime: 2019-10-24 10:27:59
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (307->27), mult. (758->66), div. (186->13), fcn. (854->7), ass. (0->41)
	t62 = sin(qJ(2));
	t56 = t62 ^ 2;
	t63 = cos(qJ(2));
	t69 = t63 ^ 2;
	t81 = t56 / t69;
	t74 = 0.1e1 + t81;
	t60 = sin(pkin(8));
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
	t61 = cos(pkin(8));
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
	t1 = [0, t35, 0, 0, 0; 0, ((-0.2e1 * t39 * t87 + (-qJD(2) * t34 - t33) * t85) * t63 + (t40 * t87 * t88 + (t35 * t40 * t76 + t86 * t88 - qJD(2) * t39 - (-t78 * qJD(2) + t75 * t43) * t46 * t84) * t36 - (t35 * t46 + (-t75 * qJD(2) + t78 * t43) * t47) * t63 * t85) * t62) * t61, 0, 0, 0; 0, (-0.1e1 / t51 - 0.2e1 / t53 / t51 ^ 2 * t82) * t61 / t60 * t77, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:59
	% EndTime: 2019-10-24 10:28:00
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (355->40), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->55)
	t97 = cos(qJ(2));
	t91 = t97 ^ 2;
	t95 = sin(qJ(2));
	t119 = 0.1e1 / t95 ^ 2 * t91;
	t110 = 0.1e1 + t119;
	t92 = sin(pkin(8));
	t128 = t110 * t92;
	t115 = qJD(2) * t97;
	t93 = cos(pkin(8));
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
	t1 = [0, t56, 0, 0, 0; 0, ((t61 * t114 + (qJD(2) * t54 + t53) * t62 * t57) * t95 + (t54 * t62 * t114 + (0.2e1 * t54 * t125 - qJD(2) * t61 - (-t56 * t78 + (t111 * qJD(2) + t116 * t60) * t79) * t123 + (t56 * t113 - (t116 * qJD(2) + t111 * t60) * t97 * t78) * t62) * t57) * t97) * t93, 0, 0, 0; 0, (0.2e1 * t109 * t97 * t124 + (t109 * t95 * qJD(2) + ((qJD(4) * t73 - 0.2e1 * t106 * t122) * t94 + (-t70 * t94 + (t69 + t127) * t96) * t74) * t97) * t65) * t93, 0, -0.2e1 * t124 - 0.2e1 * (t65 * t70 * t74 - (-t65 * t122 - t74 * t124) * t106) * t106, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:59
	% EndTime: 2019-10-24 10:28:00
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (714->44), mult. (1346->105), div. (247->12), fcn. (1511->9), ass. (0->58)
	t125 = sin(pkin(8));
	t128 = cos(qJ(2));
	t123 = t128 ^ 2;
	t127 = sin(qJ(2));
	t148 = 0.1e1 / t127 ^ 2 * t123;
	t139 = 0.1e1 + t148;
	t158 = t125 * t139;
	t145 = qJD(2) * t128;
	t147 = t125 * t128;
	t109 = atan2(-t147, t127);
	t107 = sin(t109);
	t108 = cos(t109);
	t98 = -t107 * t147 + t108 * t127;
	t95 = 0.1e1 / t98;
	t124 = qJ(4) + qJ(5);
	t115 = cos(t124);
	t114 = sin(t124);
	t126 = cos(pkin(8));
	t146 = t126 * t127;
	t142 = t114 * t146;
	t106 = t125 * t115 + t142;
	t102 = 0.1e1 / t106;
	t103 = 0.1e1 / t106 ^ 2;
	t96 = 0.1e1 / t98 ^ 2;
	t141 = t115 * t146;
	t105 = t125 * t114 - t141;
	t101 = t105 ^ 2;
	t151 = t103 * t105;
	t118 = qJD(4) + qJD(5);
	t137 = -t118 * t125 + t126 * t145;
	t92 = t137 * t114 + t118 * t141;
	t154 = t102 * t103 * t92;
	t91 = t101 * t103 + 0.1e1;
	t93 = -t137 * t115 + t118 * t142;
	t156 = (-t101 * t154 + t93 * t151) / t91 ^ 2;
	t143 = t108 * t147;
	t149 = t108 * t128;
	t150 = t107 * t127;
	t116 = t125 ^ 2;
	t113 = t116 * t148 + 0.1e1;
	t110 = 0.1e1 / t113;
	t99 = t110 * t158;
	t94 = qJD(2) * t99;
	t83 = (-t143 - t150) * t94 + (t125 * t150 + t149) * qJD(2);
	t155 = t83 * t95 * t96;
	t153 = t127 * t96;
	t152 = t125 - t99;
	t117 = t126 ^ 2;
	t88 = t117 * t123 * t96 + 0.1e1;
	t144 = 0.2e1 * (-t123 * t155 - t145 * t153) * t117 / t88 ^ 2;
	t140 = t125 * t99 - 0.1e1;
	t138 = t102 * t115 + t114 * t151;
	t89 = 0.1e1 / t91;
	t86 = 0.1e1 / t88;
	t85 = -0.2e1 * (t110 - t139 / t113 ^ 2 * t116) / t127 * t145 * t158;
	t84 = -t140 * t149 + t152 * t150;
	t80 = -0.2e1 * t156 + 0.2e1 * (t103 * t89 * t93 + (-t103 * t156 - t89 * t154) * t105) * t105;
	t1 = [0, t85, 0, 0, 0; 0, ((t95 * t144 + (qJD(2) * t84 + t83) * t96 * t86) * t127 + (t84 * t96 * t144 + (0.2e1 * t84 * t155 - qJD(2) * t95 - (-t107 * t85 + (-qJD(2) + (0.2e1 * t125 - t99) * t94) * t108) * t153 + (t85 * t143 - (t152 * qJD(2) + t140 * t94) * t128 * t107) * t96) * t86) * t128) * t126, 0, 0, 0; 0, (t138 * t89 * t127 * qJD(2) + (0.2e1 * t138 * t156 + ((t102 * t118 + 0.2e1 * t105 * t154) * t114 + (-t114 * t93 + (-t105 * t118 + t92) * t115) * t103) * t89) * t128) * t126, 0, t80, t80;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end