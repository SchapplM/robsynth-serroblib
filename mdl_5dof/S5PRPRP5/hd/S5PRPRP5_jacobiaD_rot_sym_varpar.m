% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPRP5
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
%   Wie in S5PRPRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:24
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRPRP5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:49
	% EndTime: 2019-10-24 10:24:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:49
	% EndTime: 2019-10-24 10:24:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:48
	% EndTime: 2019-10-24 10:24:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:49
	% EndTime: 2019-10-24 10:24:49
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (346->34), mult. (901->88), div. (202->12), fcn. (1017->9), ass. (0->52)
	t80 = sin(pkin(7));
	t83 = sin(qJ(2));
	t75 = t83 ^ 2;
	t84 = cos(qJ(2));
	t102 = t75 / t84 ^ 2;
	t95 = 0.1e1 + t102;
	t113 = t80 * t95;
	t82 = cos(pkin(7));
	t73 = t82 ^ 2;
	t103 = t73 * t75;
	t112 = 0.2e1 * t103;
	t111 = qJD(2) * t83;
	t101 = t80 * t83;
	t68 = atan2(-t101, -t84);
	t66 = sin(t68);
	t67 = cos(t68);
	t53 = -t66 * t101 - t67 * t84;
	t50 = 0.1e1 / t53;
	t100 = t82 * t84;
	t79 = sin(pkin(8));
	t81 = cos(pkin(8));
	t65 = t81 * t100 + t80 * t79;
	t61 = 0.1e1 / t65;
	t51 = 0.1e1 / t53 ^ 2;
	t62 = 0.1e1 / t65 ^ 2;
	t72 = t80 ^ 2;
	t71 = t72 * t102 + 0.1e1;
	t69 = 0.1e1 / t71;
	t58 = t69 * t113;
	t104 = t66 * t84;
	t93 = -t80 * t104 + t67 * t83;
	t98 = t67 * t101;
	t94 = -t98 + t104;
	t44 = t94 * t58 + t93;
	t109 = 0.2e1 * t44;
	t106 = t51 * t84;
	t54 = qJD(2) * t58;
	t43 = t93 * qJD(2) + t94 * t54;
	t107 = t43 * t50 * t51;
	t49 = t51 * t103 + 0.1e1;
	t108 = (t106 * t111 - t75 * t107) * t73 / t49 ^ 2;
	t97 = t79 * t100;
	t64 = -t80 * t81 + t97;
	t105 = t64 * t81;
	t99 = t58 - t80;
	t96 = t58 * t80 - 0.1e1;
	t63 = t61 * t62;
	t60 = t64 ^ 2;
	t57 = t60 * t62 + 0.1e1;
	t47 = 0.1e1 / t49;
	t45 = 0.2e1 * (t69 - t95 / t71 ^ 2 * t72) / t84 * t111 * t113;
	t1 = [0, t45, 0, 0, 0; 0, ((-0.2e1 * t50 * t108 + (-qJD(2) * t44 - t43) * t51 * t47) * t84 + (t51 * t108 * t109 + (t107 * t109 - qJD(2) * t50 - (t45 * t66 + (-t96 * qJD(2) + t99 * t54) * t67) * t106 + (t45 * t98 - (-t99 * qJD(2) + t96 * t54) * t83 * t66) * t51) * t47) * t83) * t82, 0, 0, 0; 0, ((-t62 * t105 + t61 * t79) / t57 ^ 2 * (t60 * t63 * t81 - t62 * t64 * t79) * t112 + (-t61 * t97 + (t63 * t105 * t112 + (t64 * t100 - 0.2e1 * t79 * t103) * t62) * t81) / t57) * qJD(2), 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:49
	% EndTime: 2019-10-24 10:24:49
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (614->45), mult. (1166->104), div. (229->12), fcn. (1315->9), ass. (0->57)
	t102 = sin(pkin(7));
	t105 = cos(qJ(2));
	t104 = sin(qJ(2));
	t98 = t104 ^ 2;
	t128 = 0.1e1 / t105 ^ 2 * t98;
	t116 = 0.1e1 + t128;
	t138 = t102 * t116;
	t137 = qJD(2) * t104;
	t103 = cos(pkin(7));
	t121 = t103 * t105;
	t96 = pkin(8) + qJ(4);
	t92 = sin(t96);
	t93 = cos(t96);
	t82 = t102 * t92 + t93 * t121;
	t122 = t102 * t104;
	t86 = atan2(-t122, -t105);
	t83 = sin(t86);
	t84 = cos(t86);
	t71 = -t84 * t105 - t83 * t122;
	t68 = 0.1e1 / t71;
	t78 = 0.1e1 / t82;
	t124 = t84 * t104;
	t125 = t105 * t83;
	t94 = t102 ^ 2;
	t89 = t94 * t128 + 0.1e1;
	t87 = 0.1e1 / t89;
	t75 = t87 * t138;
	t60 = t75 * t125 + t124 + (-t75 * t124 - t125) * t102;
	t135 = 0.2e1 * t60;
	t69 = 0.1e1 / t71 ^ 2;
	t79 = 0.1e1 / t82 ^ 2;
	t120 = qJD(2) * t105;
	t126 = t104 * t69;
	t118 = t84 * t122;
	t72 = qJD(2) * t75;
	t59 = (-t118 + t125) * t72 + (-t102 * t125 + t124) * qJD(2);
	t132 = t59 * t68 * t69;
	t95 = t103 ^ 2;
	t64 = t95 * t98 * t69 + 0.1e1;
	t134 = (t120 * t126 - t98 * t132) * t95 / t64 ^ 2;
	t81 = -t102 * t93 + t92 * t121;
	t129 = t79 * t81;
	t115 = t103 * t137;
	t123 = qJD(4) * t81;
	t74 = -t93 * t115 - t123;
	t130 = t74 * t78 * t79;
	t77 = t81 ^ 2;
	t67 = t77 * t79 + 0.1e1;
	t73 = -t82 * qJD(4) + t92 * t115;
	t133 = (-t73 * t129 - t77 * t130) / t67 ^ 2;
	t62 = 0.1e1 / t64;
	t131 = t62 * t69;
	t119 = -0.2e1 * t133;
	t114 = t93 * t129 - t78 * t92;
	t65 = 0.1e1 / t67;
	t61 = 0.2e1 * (t87 - t116 / t89 ^ 2 * t94) / t105 * t137 * t138;
	t1 = [0, t61, 0, 0, 0; 0, ((-0.2e1 * t68 * t134 + (-qJD(2) * t60 - t59) * t131) * t105 + (t69 * t134 * t135 + (t61 * t69 * t118 + t132 * t135 - qJD(2) * t68 - ((t102 * t75 - 0.1e1) * t72 + (t102 - t75) * qJD(2)) * t83 * t126) * t62 - (t61 * t83 + (qJD(2) + (-0.2e1 * t102 + t75) * t72) * t84) * t105 * t131) * t104) * t103, 0, 0, 0; 0, (t114 * t65 * t120 + (t114 * t119 + ((-qJD(4) * t78 - 0.2e1 * t81 * t130) * t93 + (-t73 * t93 + (t74 - t123) * t92) * t79) * t65) * t104) * t103, 0, t119 + 0.2e1 * (-t65 * t73 * t79 + (-t65 * t130 - t79 * t133) * t81) * t81, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:49
	% EndTime: 2019-10-24 10:24:49
	% DurationCPUTime: 0.71s
	% Computational Cost: add. (3088->84), mult. (3936->200), div. (821->15), fcn. (5024->9), ass. (0->85)
	t129 = pkin(8) + qJ(4);
	t126 = sin(t129);
	t127 = cos(t129);
	t135 = cos(pkin(7));
	t134 = sin(pkin(7));
	t137 = cos(qJ(2));
	t165 = t134 * t137;
	t113 = t126 * t165 + t135 * t127;
	t136 = sin(qJ(2));
	t163 = t136 * t126;
	t103 = atan2(-t113, t163);
	t100 = cos(t103);
	t108 = t113 ^ 2;
	t124 = 0.1e1 / t126 ^ 2;
	t132 = 0.1e1 / t136 ^ 2;
	t166 = t124 * t132;
	t104 = t108 * t166 + 0.1e1;
	t101 = 0.1e1 / t104;
	t123 = 0.1e1 / t126;
	t154 = t123 * t132 * t137;
	t147 = t113 * t154 + t134;
	t89 = t147 * t101;
	t170 = t134 - t89;
	t99 = sin(t103);
	t182 = t170 * t99 * t136 + t100 * t137;
	t164 = t135 * t137;
	t117 = t134 * t126 + t127 * t164;
	t111 = 0.1e1 / t117 ^ 2;
	t128 = t135 ^ 2;
	t130 = t136 ^ 2;
	t153 = t128 * t130 * t111;
	t107 = 0.1e1 + t153;
	t161 = qJD(2) * t137;
	t150 = t136 * t161;
	t110 = 0.1e1 / t117;
	t116 = t126 * t164 - t134 * t127;
	t162 = qJD(2) * t136;
	t151 = t135 * t162;
	t98 = -t116 * qJD(4) - t127 * t151;
	t175 = t110 * t111 * t98;
	t156 = t130 * t175;
	t181 = 0.1e1 / t107 ^ 2 * (t111 * t150 - t156) * t128;
	t131 = 0.1e1 / t136;
	t115 = -t135 * t126 + t127 * t165;
	t167 = t124 * t127;
	t146 = t113 * t167 - t115 * t123;
	t180 = t131 * t146;
	t125 = t123 * t124;
	t133 = t131 / t130;
	t155 = t113 * t166;
	t159 = qJD(4) * t127;
	t148 = t137 * t159;
	t152 = t134 * t162;
	t95 = -t134 * t148 + (qJD(4) * t135 + t152) * t126;
	t179 = -0.2e1 / t104 ^ 2 * (-t95 * t155 + (-t124 * t133 * t161 - t125 * t132 * t159) * t108);
	t171 = t99 * t113;
	t94 = t100 * t163 - t171;
	t90 = 0.1e1 / t94;
	t91 = 0.1e1 / t94 ^ 2;
	t178 = 0.2e1 * t116;
	t149 = t126 * t161;
	t145 = t136 * t159 + t149;
	t83 = (t95 * t131 * t123 + t145 * t155) * t101;
	t174 = t113 * t83;
	t80 = (-t83 * t163 + t95) * t99 + (t145 - t174) * t100;
	t92 = t90 * t91;
	t177 = t80 * t92;
	t160 = qJD(4) * t126;
	t97 = t126 * t151 - t134 * t160 - t135 * t148;
	t176 = t97 * t91;
	t173 = t116 * t91;
	t169 = t100 * t113;
	t109 = t116 ^ 2;
	t88 = t109 * t91 + 0.1e1;
	t158 = 0.2e1 * (-t109 * t177 - t97 * t173) / t88 ^ 2;
	t157 = t136 * t178;
	t105 = 0.1e1 / t107;
	t96 = t113 * qJD(4) + t127 * t152;
	t86 = 0.1e1 / t88;
	t85 = t101 * t180;
	t82 = t182 * t126 - t89 * t169;
	t81 = (-t85 * t163 - t115) * t99 + (-t113 * t85 + t127 * t136) * t100;
	t79 = t147 * t179 + (-t95 * t154 + (-t148 * t166 + (-0.2e1 * t133 * t137 ^ 2 - t131) * t123 * qJD(2)) * t113) * t101;
	t77 = t179 * t180 + (-t146 * t132 * t161 + (-t95 * t167 + t123 * t96 + (t115 * t167 + (-0.2e1 * t125 * t127 ^ 2 - t123) * t113) * qJD(4)) * t131) * t101;
	t1 = [0, t79, 0, t77, 0; 0, t82 * t158 * t173 + (t82 * t176 + (-(-t79 * t169 + (t100 * t95 + t171 * t83) * t89) * t91 + 0.2e1 * t82 * t177) * t116 + (-t135 * t136 * t90 - t182 * t173) * t159) * t86 + (-((t170 * qJD(2) - t83) * t99 * t137 + (-t79 * t99 + (t170 * t83 - qJD(2)) * t100) * t136) * t86 * t173 + (-t90 * t86 * t161 + (t91 * t86 * t80 + t90 * t158) * t136) * t135) * t126, 0, (-t117 * t90 + t81 * t173) * t158 + (t81 * t176 + t98 * t90 + (t81 * t92 * t178 - t117 * t91) * t80 + (-(t96 + (-t149 + t174) * t85 + (-t126 * t77 + (-qJD(4) * t85 - t83) * t127) * t136) * t99 - (t127 * t161 - t113 * t77 - t115 * t83 + t85 * t95 + (-t83 * t85 - qJD(4)) * t163) * t100) * t173) * t86, 0; 0, 0.2e1 * (t110 * t164 + t127 * t153) * t181 + (0.2e1 * t127 * t128 * t156 + t110 * t151 + (t98 * t164 + (-0.2e1 * t127 * t150 + t130 * t160) * t128) * t111) * t105, 0, (t111 * t157 * t181 + (t157 * t175 + (-t116 * t161 + t136 * t97) * t111) * t105) * t135, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end