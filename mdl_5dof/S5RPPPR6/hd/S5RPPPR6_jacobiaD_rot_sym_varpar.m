% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPPR6
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
%   Wie in S5RPPPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPPPR6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:54:08
	% EndTime: 2019-12-29 15:54:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:54:08
	% EndTime: 2019-12-29 15:54:08
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:54:08
	% EndTime: 2019-12-29 15:54:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:54:08
	% EndTime: 2019-12-29 15:54:08
	% DurationCPUTime: 0.48s
	% Computational Cost: add. (213->24), mult. (468->72), div. (98->14), fcn. (596->7), ass. (0->34)
	t70 = cos(qJ(1));
	t89 = 0.2e1 * t70;
	t67 = sin(pkin(7));
	t56 = t67 ^ 2;
	t68 = cos(pkin(7));
	t59 = 0.1e1 / t68 ^ 2;
	t69 = sin(qJ(1));
	t61 = t69 ^ 2;
	t53 = t56 * t59 * t61 + 0.1e1;
	t49 = 0.1e1 / t53;
	t87 = (t49 - 0.1e1) * t67;
	t81 = t69 * t67;
	t48 = atan2(-t81, -t68);
	t46 = sin(t48);
	t47 = cos(t48);
	t44 = -t46 * t81 - t47 * t68;
	t41 = 0.1e1 / t44;
	t58 = 0.1e1 / t68;
	t42 = 0.1e1 / t44 ^ 2;
	t85 = t42 * t70;
	t66 = t70 ^ 2;
	t84 = 0.1e1 / t53 ^ 2 * t66;
	t83 = t56 * t58;
	t82 = 0.1e1 / t69 ^ 2 * t66;
	t80 = qJD(1) * t69;
	t57 = t68 ^ 2;
	t79 = t58 / t57 * t84;
	t37 = (-t47 * t49 * t69 * t83 + t46 * t87) * t70;
	t55 = t67 * t56;
	t54 = t57 * t82 + 0.1e1;
	t43 = t41 * t42;
	t40 = t42 * t56 * t66 + 0.1e1;
	t36 = qJD(1) * t37;
	t1 = [(-t49 * t58 * t67 - 0.2e1 * t55 * t79) * t80, 0, 0, 0, 0; (0.2e1 * (t37 * t85 + t41 * t69) / t40 ^ 2 * (-t36 * t43 * t66 - t80 * t85) * t56 + ((t37 * t43 * t89 + t42 * t69) * t36 + (-t70 * t41 + ((t37 + (t55 * t59 * t84 + t87) * t70 * t46) * t69 - (0.2e1 * t61 * t56 ^ 2 * t79 + (t84 + (t61 - 0.2e1 * t66) * t49) * t83) * t70 * t47) * t42) * qJD(1)) / t40) * t67, 0, 0, 0, 0; (0.1e1 / t54 - (0.1e1 + t82) / t54 ^ 2 * t57) * (0.1e1 + 0.1e1 / t61 * t66) / t69 * qJD(1) * t68 * t89, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:54:08
	% EndTime: 2019-12-29 15:54:09
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (201->31), mult. (614->96), div. (108->12), fcn. (792->9), ass. (0->49)
	t85 = sin(pkin(7));
	t78 = 0.1e1 / t85 ^ 2;
	t87 = cos(pkin(7));
	t81 = t87 ^ 2;
	t88 = sin(qJ(1));
	t82 = t88 ^ 2;
	t76 = t78 * t81 * t82 + 0.1e1;
	t74 = 0.1e1 / t76;
	t112 = (t74 - 0.1e1) * t87;
	t101 = t88 * t87;
	t73 = atan2(-t101, t85);
	t71 = sin(t73);
	t72 = cos(t73);
	t57 = -t71 * t101 + t72 * t85;
	t54 = 0.1e1 / t57;
	t86 = cos(pkin(8));
	t102 = t88 * t86;
	t84 = sin(pkin(8));
	t89 = cos(qJ(1));
	t105 = t84 * t89;
	t68 = t85 * t105 + t102;
	t64 = 0.1e1 / t68;
	t77 = 0.1e1 / t85;
	t55 = 0.1e1 / t57 ^ 2;
	t65 = 0.1e1 / t68 ^ 2;
	t110 = t55 * t89;
	t103 = t88 * t84;
	t104 = t86 * t89;
	t98 = t85 * t104 - t103;
	t109 = t65 * t98;
	t70 = -t85 * t103 + t104;
	t108 = t98 * t70;
	t83 = t89 ^ 2;
	t107 = 0.1e1 / t76 ^ 2 * t83;
	t106 = t77 * t81;
	t100 = qJD(1) * t88;
	t99 = t87 * t81 * t107;
	t69 = t85 * t102 + t105;
	t50 = (t72 * t74 * t88 * t106 + t71 * t112) * t89;
	t79 = t77 * t78;
	t66 = t64 * t65;
	t63 = t98 ^ 2;
	t62 = t70 * qJD(1);
	t61 = t69 * qJD(1);
	t60 = t63 * t65 + 0.1e1;
	t56 = t54 * t55;
	t53 = t55 * t81 * t83 + 0.1e1;
	t49 = qJD(1) * t50;
	t1 = [(t74 * t77 * t87 + 0.2e1 * t79 * t99) * t100, 0, 0, 0, 0; (0.2e1 * (t50 * t110 + t54 * t88) / t53 ^ 2 * (-t49 * t56 * t83 - t100 * t110) * t81 + ((0.2e1 * t50 * t56 * t89 + t55 * t88) * t49 + (-t89 * t54 + ((t50 + (t78 * t99 + t112) * t89 * t71) * t88 - (-0.2e1 * t79 * t82 * t81 ^ 2 * t107 + (-t107 + (-t82 + 0.2e1 * t83) * t74) * t106) * t89 * t72) * t55) * qJD(1)) / t53) * t87, 0, 0, 0, 0; 0.2e1 * (-t65 * t108 - t64 * t69) / t60 ^ 2 * (-t62 * t63 * t66 - t61 * t109) + (-t70 * t61 * t65 + (-0.2e1 * t66 * t108 - t69 * t65) * t62 + (-t68 * t109 + t98 * t64) * qJD(1)) / t60, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:54:09
	% EndTime: 2019-12-29 15:54:09
	% DurationCPUTime: 0.71s
	% Computational Cost: add. (496->44), mult. (1652->112), div. (172->13), fcn. (2153->11), ass. (0->67)
	t146 = cos(pkin(8));
	t147 = cos(pkin(7));
	t177 = t147 * t146;
	t148 = sin(qJ(5));
	t150 = cos(qJ(5));
	t145 = sin(pkin(8));
	t149 = sin(qJ(1));
	t186 = sin(pkin(7));
	t189 = cos(qJ(1));
	t163 = t189 * t186;
	t159 = -t145 * t163 - t149 * t146;
	t167 = t147 * t189;
	t161 = t148 * t159 + t150 * t167;
	t113 = t161 ^ 2;
	t118 = t148 * t167 - t150 * t159;
	t115 = 0.1e1 / t118 ^ 2;
	t112 = t113 * t115 + 0.1e1;
	t110 = 0.1e1 / t112;
	t165 = t149 * t186;
	t135 = -t145 * t165 + t189 * t146;
	t130 = t135 * qJD(1);
	t175 = qJD(1) * t149;
	t166 = t147 * t175;
	t109 = t161 * qJD(5) + t130 * t150 - t148 * t166;
	t114 = 0.1e1 / t118;
	t184 = t109 * t114 * t115;
	t108 = t118 * qJD(5) + t130 * t148 + t150 * t166;
	t182 = t115 * t161;
	t187 = 0.1e1 / t112 ^ 2 * (-t108 * t182 - t113 * t184);
	t192 = t161 * (t110 * t184 + t115 * t187);
	t179 = 0.1e1 / t177;
	t178 = 0.1e1 / t146 ^ 2 / t147 ^ 2;
	t134 = t189 * t145 + t146 * t165;
	t125 = atan2(t134, t177);
	t121 = sin(t125);
	t122 = cos(t125);
	t107 = t121 * t134 + t122 * t177;
	t104 = 0.1e1 / t107;
	t190 = t134 ^ 2;
	t105 = 0.1e1 / t107 ^ 2;
	t162 = t146 * t163;
	t127 = -qJD(1) * t162 + t145 * t175;
	t126 = t190 * t178 + 0.1e1;
	t123 = 0.1e1 / t126;
	t160 = -t121 + (-t122 * t134 * t179 + t121) * t123;
	t99 = t160 * t127;
	t188 = t104 * t105 * t99;
	t132 = t145 * t149 - t162;
	t185 = t105 * t132;
	t183 = t110 * t115;
	t124 = 0.1e1 / t126 ^ 2;
	t181 = t124 * t134;
	t180 = t127 * t132;
	t176 = t147 * t149;
	t174 = -0.2e1 * t187;
	t172 = -0.2e1 * t178 * t179;
	t170 = t108 * t183;
	t168 = t123 * t179;
	t164 = qJD(1) * t167;
	t119 = t135 * t148 + t150 * t176;
	t120 = t135 * t150 - t148 * t176;
	t131 = t132 ^ 2;
	t129 = t134 * qJD(1);
	t128 = t159 * qJD(1);
	t103 = t105 * t131 + 0.1e1;
	t100 = t160 * t132;
	t1 = [t172 * t180 * t181 - t129 * t168, 0, 0, 0, 0; 0.2e1 * (t100 * t185 - t104 * t134) * (t129 * t185 - t131 * t188) / t103 ^ 2 + (-t127 * t104 + (-t100 * t129 - t134 * t99) * t105 + (0.2e1 * t100 * t188 + (-t160 * t129 - (t121 * t178 * t181 + (0.2e1 * t168 + (t190 * t172 - t179) * t124) * t122) * t180) * t105) * t132) / t103, 0, 0, 0, 0; -t120 * t170 - 0.2e1 * t120 * t192 + (-t109 * t183 + t114 * t174) * t119 + ((t120 * qJD(5) + t128 * t148 + t150 * t164) * t114 + (-t119 * qJD(5) + t128 * t150 - t148 * t164) * t182) * t110, 0, 0, 0, t174 - 0.2e1 * (t170 + t192) * t161;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end