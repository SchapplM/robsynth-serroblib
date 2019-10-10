% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6PRRRRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:09
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRP5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:09:35
	% EndTime: 2019-10-09 23:09:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:09:35
	% EndTime: 2019-10-09 23:09:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:09:35
	% EndTime: 2019-10-09 23:09:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:09:35
	% EndTime: 2019-10-09 23:09:35
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (179->22), mult. (525->56), div. (35->9), fcn. (736->13), ass. (0->37)
	t45 = sin(pkin(12));
	t48 = cos(pkin(12));
	t54 = cos(qJ(2));
	t50 = cos(pkin(6));
	t52 = sin(qJ(2));
	t58 = t50 * t52;
	t43 = -t45 * t58 + t48 * t54;
	t51 = sin(qJ(3));
	t63 = t43 * t51;
	t53 = cos(qJ(3));
	t62 = t43 * t53;
	t46 = sin(pkin(7));
	t47 = sin(pkin(6));
	t61 = t46 * t47;
	t49 = cos(pkin(7));
	t60 = t47 * t49;
	t59 = t47 * t52;
	t57 = t50 * t54;
	t42 = -t45 * t57 - t48 * t52;
	t55 = t42 * t49 + t45 * t61;
	t32 = t55 * t51 + t62;
	t30 = 0.1e1 / t32 ^ 2;
	t31 = -t55 * t53 + t63;
	t56 = t31 ^ 2 * t30 + 0.1e1;
	t41 = -t45 * t54 - t48 * t58;
	t40 = t50 * t49 - t54 * t61;
	t39 = 0.1e1 / t40 ^ 2;
	t38 = -t42 * t46 + t45 * t60;
	t37 = (-t45 * t52 + t48 * t57) * t46 + t48 * t60;
	t36 = atan2(t37, t40);
	t34 = cos(t36);
	t33 = sin(t36);
	t29 = 0.1e1 / t56;
	t28 = t33 * t37 + t34 * t40;
	t27 = 0.1e1 / t28 ^ 2;
	t25 = (t41 / t40 - t37 * t39 * t59) * t46 / (t37 ^ 2 * t39 + 0.1e1);
	t1 = [0, t25, 0, 0, 0, 0; 0, (t43 * t46 / t28 - ((t33 * t41 + t34 * t59) * t46 + (-t33 * t40 + t34 * t37) * t25) * t38 * t27) / (t38 ^ 2 * t27 + 0.1e1), 0, 0, 0, 0; 0, ((t42 * t51 + t49 * t62) / t32 - (t42 * t53 - t49 * t63) * t31 * t30) * t29, t56 * t29, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:09:35
	% EndTime: 2019-10-09 23:09:35
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (607->40), mult. (1825->97), div. (65->9), fcn. (2498->15), ass. (0->60)
	t87 = sin(pkin(7));
	t88 = sin(pkin(6));
	t109 = t88 * t87;
	t89 = cos(pkin(12));
	t90 = cos(pkin(7));
	t91 = cos(pkin(6));
	t97 = cos(qJ(2));
	t106 = t91 * t97;
	t86 = sin(pkin(12));
	t94 = sin(qJ(2));
	t99 = t89 * t106 - t86 * t94;
	t116 = -t89 * t109 + t99 * t90;
	t107 = t91 * t94;
	t81 = t89 * t107 + t86 * t97;
	t93 = sin(qJ(3));
	t96 = cos(qJ(3));
	t66 = -t116 * t96 + t81 * t93;
	t108 = t90 * t96;
	t110 = t87 * t91;
	t75 = -t96 * t110 + (-t97 * t108 + t93 * t94) * t88;
	t65 = atan2(-t66, t75);
	t62 = sin(t65);
	t63 = cos(t65);
	t56 = -t62 * t66 + t63 * t75;
	t55 = 0.1e1 / t56 ^ 2;
	t103 = t86 * t109;
	t83 = -t86 * t107 + t89 * t97;
	t111 = t83 * t93;
	t82 = -t86 * t106 - t89 * t94;
	t69 = -t96 * t103 - t82 * t108 + t111;
	t115 = t55 * t69;
	t70 = t83 * t96 + (t82 * t90 + t103) * t93;
	t77 = t86 * t88 * t90 - t82 * t87;
	t92 = sin(qJ(4));
	t95 = cos(qJ(4));
	t61 = t70 * t95 + t77 * t92;
	t59 = 0.1e1 / t61 ^ 2;
	t60 = t70 * t92 - t77 * t95;
	t114 = t59 * t60;
	t74 = 0.1e1 / t75 ^ 2;
	t113 = t66 * t74;
	t112 = t83 * t87;
	t105 = t93 * t97;
	t104 = t94 * t96;
	t101 = t60 ^ 2 * t59 + 0.1e1;
	t100 = -t62 * t75 - t63 * t66;
	t80 = (t90 * t104 + t105) * t88;
	t76 = t93 * t110 + (t90 * t105 + t104) * t88;
	t73 = 0.1e1 / t75;
	t72 = -t90 * t111 + t82 * t96;
	t71 = t81 * t108 + t99 * t93;
	t68 = t116 * t93 + t81 * t96;
	t64 = 0.1e1 / (t66 ^ 2 * t74 + 0.1e1);
	t58 = 0.1e1 / t61;
	t57 = 0.1e1 / t101;
	t54 = 0.1e1 / t56;
	t53 = 0.1e1 / (t69 ^ 2 * t55 + 0.1e1);
	t52 = (t80 * t113 - t71 * t73) * t64;
	t51 = (t76 * t113 - t68 * t73) * t64;
	t1 = [0, t52, t51, 0, 0, 0; 0, ((t83 * t108 + t82 * t93) * t54 - (t100 * t52 - t62 * t71 + t63 * t80) * t115) * t53, (t70 * t54 - (t100 * t51 - t62 * t68 + t63 * t76) * t115) * t53, 0, 0, 0; 0, ((-t95 * t112 + t72 * t92) * t58 - (t92 * t112 + t72 * t95) * t114) * t57, (t95 * t114 - t92 * t58) * t69 * t57, t101 * t57, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:09:35
	% EndTime: 2019-10-09 23:09:36
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (1524->59), mult. (4378->143), div. (95->9), fcn. (6002->17), ass. (0->75)
	t122 = sin(pkin(6));
	t128 = sin(qJ(3));
	t129 = sin(qJ(2));
	t132 = cos(qJ(3));
	t133 = cos(qJ(2));
	t124 = cos(pkin(7));
	t146 = t124 * t128;
	t121 = sin(pkin(7));
	t125 = cos(pkin(6));
	t150 = t121 * t125;
	t112 = t128 * t150 + (t129 * t132 + t133 * t146) * t122;
	t148 = t122 * t121;
	t114 = t125 * t124 - t133 * t148;
	t127 = sin(qJ(4));
	t131 = cos(qJ(4));
	t104 = t112 * t127 - t114 * t131;
	t120 = sin(pkin(12));
	t123 = cos(pkin(12));
	t144 = t125 * t129;
	t115 = t120 * t133 + t123 * t144;
	t143 = t125 * t133;
	t136 = -t120 * t129 + t123 * t143;
	t134 = -t123 * t148 + t136 * t124;
	t101 = t115 * t132 + t134 * t128;
	t147 = t122 * t124;
	t135 = -t136 * t121 - t123 * t147;
	t91 = t101 * t127 - t135 * t131;
	t90 = atan2(-t91, t104);
	t87 = sin(t90);
	t88 = cos(t90);
	t81 = t104 * t88 - t87 * t91;
	t80 = 0.1e1 / t81 ^ 2;
	t116 = -t120 * t143 - t123 * t129;
	t117 = -t120 * t144 + t123 * t133;
	t139 = t120 * t148;
	t103 = t117 * t132 + (t116 * t124 + t139) * t128;
	t137 = -t116 * t121 + t120 * t147;
	t94 = t103 * t127 - t137 * t131;
	t154 = t80 * t94;
	t145 = t124 * t132;
	t102 = -t116 * t145 + t117 * t128 - t132 * t139;
	t126 = sin(qJ(5));
	t130 = cos(qJ(5));
	t95 = t103 * t131 + t137 * t127;
	t86 = t102 * t126 + t130 * t95;
	t84 = 0.1e1 / t86 ^ 2;
	t85 = -t102 * t130 + t126 * t95;
	t153 = t84 * t85;
	t99 = 0.1e1 / t104 ^ 2;
	t152 = t91 * t99;
	t151 = t102 * t131;
	t149 = t121 * t131;
	t142 = t128 * t129;
	t141 = t132 * t133;
	t140 = t84 * t85 ^ 2 + 0.1e1;
	t138 = -t104 * t87 - t88 * t91;
	t111 = t132 * t150 + (t124 * t141 - t142) * t122;
	t108 = ((-t124 * t142 + t141) * t127 - t129 * t149) * t122;
	t107 = t116 * t132 - t117 * t146;
	t106 = t116 * t128 + t117 * t145;
	t105 = t112 * t131 + t114 * t127;
	t100 = -t115 * t128 + t134 * t132;
	t98 = 0.1e1 / t104;
	t97 = t117 * t121 * t127 + t107 * t131;
	t96 = (-t115 * t146 + t136 * t132) * t127 - t115 * t149;
	t93 = t101 * t131 + t135 * t127;
	t89 = 0.1e1 / (t91 ^ 2 * t99 + 0.1e1);
	t83 = 0.1e1 / t86;
	t82 = 0.1e1 / t140;
	t79 = 0.1e1 / t81;
	t78 = 0.1e1 / (t80 * t94 ^ 2 + 0.1e1);
	t77 = (-t100 * t98 + t111 * t152) * t89 * t127;
	t76 = (t108 * t152 - t96 * t98) * t89;
	t75 = (t105 * t152 - t93 * t98) * t89;
	t1 = [0, t76, t77, t75, 0, 0; 0, ((t107 * t127 - t117 * t149) * t79 - (t108 * t88 + t138 * t76 - t87 * t96) * t154) * t78, (-t102 * t127 * t79 - (t138 * t77 + (-t100 * t87 + t111 * t88) * t127) * t154) * t78, (t95 * t79 - (t105 * t88 + t138 * t75 - t87 * t93) * t154) * t78, 0, 0; 0, ((-t106 * t130 + t97 * t126) * t83 - (t106 * t126 + t97 * t130) * t153) * t82, ((-t103 * t130 - t126 * t151) * t83 - (t103 * t126 - t130 * t151) * t153) * t82, (-t126 * t83 + t130 * t153) * t94 * t82, t140 * t82, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:09:35
	% EndTime: 2019-10-09 23:09:36
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (1524->59), mult. (4378->143), div. (95->9), fcn. (6002->17), ass. (0->75)
	t129 = sin(pkin(12));
	t132 = cos(pkin(12));
	t138 = sin(qJ(2));
	t134 = cos(pkin(6));
	t142 = cos(qJ(2));
	t152 = t134 * t142;
	t125 = -t129 * t152 - t132 * t138;
	t153 = t134 * t138;
	t126 = -t129 * t153 + t132 * t142;
	t133 = cos(pkin(7));
	t137 = sin(qJ(3));
	t141 = cos(qJ(3));
	t130 = sin(pkin(7));
	t131 = sin(pkin(6));
	t157 = t131 * t130;
	t148 = t129 * t157;
	t112 = t126 * t141 + (t125 * t133 + t148) * t137;
	t136 = sin(qJ(4));
	t140 = cos(qJ(4));
	t156 = t131 * t133;
	t146 = -t125 * t130 + t129 * t156;
	t104 = t112 * t140 + t146 * t136;
	t154 = t133 * t141;
	t111 = -t125 * t154 + t126 * t137 - t141 * t148;
	t135 = sin(qJ(5));
	t139 = cos(qJ(5));
	t95 = t104 * t139 + t111 * t135;
	t93 = 0.1e1 / t95 ^ 2;
	t94 = t104 * t135 - t111 * t139;
	t163 = t93 * t94;
	t103 = t112 * t136 - t146 * t140;
	t124 = t129 * t142 + t132 * t153;
	t145 = -t129 * t138 + t132 * t152;
	t143 = -t132 * t157 + t145 * t133;
	t110 = t124 * t141 + t143 * t137;
	t144 = -t145 * t130 - t132 * t156;
	t100 = t110 * t136 - t144 * t140;
	t155 = t133 * t137;
	t159 = t130 * t134;
	t121 = t137 * t159 + (t138 * t141 + t142 * t155) * t131;
	t123 = t134 * t133 - t142 * t157;
	t113 = t121 * t136 - t123 * t140;
	t99 = atan2(-t100, t113);
	t96 = sin(t99);
	t97 = cos(t99);
	t90 = -t96 * t100 + t97 * t113;
	t89 = 0.1e1 / t90 ^ 2;
	t162 = t103 * t89;
	t108 = 0.1e1 / t113 ^ 2;
	t161 = t100 * t108;
	t160 = t111 * t140;
	t158 = t130 * t140;
	t151 = t137 * t138;
	t150 = t141 * t142;
	t149 = t94 ^ 2 * t93 + 0.1e1;
	t147 = -t100 * t97 - t113 * t96;
	t120 = t141 * t159 + (t133 * t150 - t151) * t131;
	t117 = ((-t133 * t151 + t150) * t136 - t138 * t158) * t131;
	t116 = t125 * t141 - t126 * t155;
	t115 = t125 * t137 + t126 * t154;
	t114 = t121 * t140 + t123 * t136;
	t109 = -t124 * t137 + t143 * t141;
	t107 = 0.1e1 / t113;
	t106 = t126 * t130 * t136 + t116 * t140;
	t105 = (-t124 * t155 + t145 * t141) * t136 - t124 * t158;
	t102 = t110 * t140 + t144 * t136;
	t98 = 0.1e1 / (t100 ^ 2 * t108 + 0.1e1);
	t92 = 0.1e1 / t95;
	t91 = 0.1e1 / t149;
	t88 = 0.1e1 / t90;
	t87 = 0.1e1 / (t103 ^ 2 * t89 + 0.1e1);
	t86 = (-t107 * t109 + t120 * t161) * t98 * t136;
	t85 = (-t105 * t107 + t117 * t161) * t98;
	t84 = (-t102 * t107 + t114 * t161) * t98;
	t1 = [0, t85, t86, t84, 0, 0; 0, ((t116 * t136 - t126 * t158) * t88 - (-t96 * t105 + t97 * t117 + t147 * t85) * t162) * t87, (-t111 * t136 * t88 - (t147 * t86 + (-t109 * t96 + t120 * t97) * t136) * t162) * t87, (t104 * t88 - (-t96 * t102 + t97 * t114 + t147 * t84) * t162) * t87, 0, 0; 0, ((t106 * t135 - t115 * t139) * t92 - (t106 * t139 + t115 * t135) * t163) * t91, ((-t112 * t139 - t135 * t160) * t92 - (t112 * t135 - t139 * t160) * t163) * t91, (-t135 * t92 + t139 * t163) * t91 * t103, t149 * t91, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end