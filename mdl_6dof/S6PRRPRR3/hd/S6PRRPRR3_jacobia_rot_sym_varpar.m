% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR3
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
%   Wie in S6PRRPRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:29
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.17s
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
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (239->26), mult. (687->64), div. (36->9), fcn. (960->15), ass. (0->42)
	t66 = sin(pkin(12));
	t67 = sin(pkin(7));
	t68 = sin(pkin(6));
	t82 = t67 * t68;
	t84 = t66 * t82;
	t71 = cos(pkin(7));
	t65 = sin(pkin(13));
	t69 = cos(pkin(13));
	t73 = sin(qJ(3));
	t75 = cos(qJ(3));
	t77 = t75 * t65 + t73 * t69;
	t55 = t77 * t71;
	t70 = cos(pkin(12));
	t74 = sin(qJ(2));
	t72 = cos(pkin(6));
	t76 = cos(qJ(2));
	t78 = t72 * t76;
	t59 = -t66 * t78 - t70 * t74;
	t79 = t72 * t74;
	t60 = -t66 * t79 + t70 * t76;
	t61 = t73 * t65 - t75 * t69;
	t46 = t59 * t55 - t60 * t61 + t77 * t84;
	t43 = 0.1e1 / t46 ^ 2;
	t54 = t61 * t71;
	t44 = -t59 * t54 - t60 * t77 - t61 * t84;
	t83 = t44 ^ 2 * t43;
	t81 = t68 * t71;
	t80 = t68 * t74;
	t58 = -t66 * t76 - t70 * t79;
	t57 = t72 * t71 - t76 * t82;
	t56 = 0.1e1 / t57 ^ 2;
	t52 = -t59 * t67 + t66 * t81;
	t51 = (-t66 * t74 + t70 * t78) * t67 + t70 * t81;
	t50 = atan2(t51, t57);
	t48 = cos(t50);
	t47 = sin(t50);
	t42 = 0.1e1 / t46;
	t41 = t47 * t51 + t48 * t57;
	t40 = 0.1e1 / t41 ^ 2;
	t38 = 0.1e1 / (0.1e1 + t83);
	t37 = (t58 / t57 - t51 * t56 * t80) * t67 / (t51 ^ 2 * t56 + 0.1e1);
	t1 = [0, t37, 0, 0, 0, 0; 0, (t60 * t67 / t41 - ((t47 * t58 + t48 * t80) * t67 + (-t47 * t57 + t48 * t51) * t37) * t52 * t40) / (t52 ^ 2 * t40 + 0.1e1), 0, 0, 0, 0; 0, ((-t60 * t54 + t59 * t77) * t42 + (-t60 * t55 - t59 * t61) * t44 * t43) * t38, (t42 * t46 + t83) * t38, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (1073->42), mult. (3064->103), div. (65->9), fcn. (4203->17), ass. (0->62)
	t106 = cos(pkin(12));
	t104 = sin(pkin(6));
	t103 = sin(pkin(7));
	t101 = sin(pkin(13));
	t105 = cos(pkin(13));
	t110 = sin(qJ(3));
	t113 = cos(qJ(3));
	t97 = t110 * t101 - t113 * t105;
	t116 = t97 * t103;
	t115 = t104 * t116;
	t117 = t113 * t101 + t110 * t105;
	t107 = cos(pkin(7));
	t91 = t97 * t107;
	t102 = sin(pkin(12));
	t111 = sin(qJ(2));
	t108 = cos(pkin(6));
	t114 = cos(qJ(2));
	t120 = t108 * t114;
	t93 = -t102 * t111 + t106 * t120;
	t121 = t108 * t111;
	t94 = t102 * t114 + t106 * t121;
	t77 = t106 * t115 - t117 * t94 - t93 * t91;
	t84 = (t111 * t117 + t114 * t91) * t104 + t108 * t116;
	t71 = atan2(t77, t84);
	t68 = sin(t71);
	t69 = cos(t71);
	t66 = t68 * t77 + t69 * t84;
	t65 = 0.1e1 / t66 ^ 2;
	t95 = -t102 * t120 - t106 * t111;
	t96 = -t102 * t121 + t106 * t114;
	t78 = -t102 * t115 - t117 * t96 - t95 * t91;
	t126 = t65 * t78;
	t109 = sin(qJ(5));
	t112 = cos(qJ(5));
	t122 = t102 * t104;
	t90 = t117 * t103;
	t92 = t117 * t107;
	t80 = t90 * t122 + t95 * t92 - t96 * t97;
	t88 = -t95 * t103 + t107 * t122;
	t75 = t88 * t109 + t80 * t112;
	t73 = 0.1e1 / t75 ^ 2;
	t74 = t80 * t109 - t88 * t112;
	t125 = t73 * t74;
	t82 = 0.1e1 / t84 ^ 2;
	t124 = t77 * t82;
	t123 = t103 * t96;
	t119 = t74 ^ 2 * t73 + 0.1e1;
	t118 = -t68 * t84 + t69 * t77;
	t87 = (-t111 * t91 + t114 * t117) * t104;
	t86 = -t96 * t92 - t95 * t97;
	t85 = -t117 * t93 + t94 * t91;
	t83 = t108 * t90 + (-t111 * t97 + t114 * t92) * t104;
	t81 = 0.1e1 / t84;
	t76 = t106 * t104 * t90 - t93 * t92 + t94 * t97;
	t72 = 0.1e1 / t75;
	t70 = 0.1e1 / (t77 ^ 2 * t82 + 0.1e1);
	t67 = 0.1e1 / t119;
	t64 = 0.1e1 / t66;
	t63 = 0.1e1 / (t78 ^ 2 * t65 + 0.1e1);
	t62 = (-t87 * t124 + t81 * t85) * t70;
	t61 = (-t83 * t124 + t76 * t81) * t70;
	t1 = [0, t62, t61, 0, 0, 0; 0, ((t117 * t95 - t96 * t91) * t64 + (t118 * t62 + t68 * t85 + t69 * t87) * t126) * t63, (t80 * t64 + (t118 * t61 + t68 * t76 + t69 * t83) * t126) * t63, 0, 0, 0; 0, ((t86 * t109 - t112 * t123) * t72 - (t109 * t123 + t86 * t112) * t125) * t67, (t109 * t72 - t112 * t125) * t78 * t67, 0, t119 * t67, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:27
	% DurationCPUTime: 0.68s
	% Computational Cost: add. (2312->62), mult. (6475->149), div. (95->9), fcn. (8887->19), ass. (0->77)
	t165 = cos(qJ(3));
	t141 = sin(qJ(5));
	t145 = cos(qJ(5));
	t133 = sin(pkin(13));
	t142 = sin(qJ(3));
	t161 = cos(pkin(13));
	t130 = -t133 * t165 - t142 * t161;
	t135 = sin(pkin(7));
	t122 = t130 * t135;
	t138 = cos(pkin(7));
	t124 = t130 * t138;
	t134 = sin(pkin(12));
	t137 = cos(pkin(12));
	t143 = sin(qJ(2));
	t139 = cos(pkin(6));
	t146 = cos(qJ(2));
	t154 = t139 * t146;
	t127 = -t134 * t154 - t137 * t143;
	t155 = t139 * t143;
	t128 = -t134 * t155 + t137 * t146;
	t149 = -t142 * t133 + t161 * t165;
	t136 = sin(pkin(6));
	t159 = t134 * t136;
	t148 = -t122 * t159 - t127 * t124 + t128 * t149;
	t156 = t136 * t138;
	t151 = -t127 * t135 + t134 * t156;
	t102 = t141 * t151 + t145 * t148;
	t121 = t149 * t135;
	t123 = t149 * t138;
	t110 = t121 * t159 + t127 * t123 + t128 * t130;
	t140 = sin(qJ(6));
	t144 = cos(qJ(6));
	t93 = t102 * t144 - t110 * t140;
	t91 = 0.1e1 / t93 ^ 2;
	t92 = t102 * t140 + t110 * t144;
	t164 = t91 * t92;
	t101 = t141 * t148 - t145 * t151;
	t115 = -t139 * t122 + (-t124 * t146 + t143 * t149) * t136;
	t125 = -t136 * t146 * t135 + t139 * t138;
	t112 = t115 * t141 - t125 * t145;
	t126 = t134 * t146 + t137 * t155;
	t150 = -t134 * t143 + t137 * t154;
	t157 = t136 * t137;
	t108 = t122 * t157 - t124 * t150 + t126 * t149;
	t147 = -t135 * t150 - t137 * t156;
	t98 = t108 * t141 - t145 * t147;
	t97 = atan2(-t98, t112);
	t94 = sin(t97);
	t95 = cos(t97);
	t88 = t112 * t95 - t94 * t98;
	t87 = 0.1e1 / t88 ^ 2;
	t163 = t101 * t87;
	t106 = 0.1e1 / t112 ^ 2;
	t162 = t106 * t98;
	t160 = t110 * t145;
	t158 = t135 * t145;
	t153 = t91 * t92 ^ 2 + 0.1e1;
	t152 = -t112 * t94 - t95 * t98;
	t118 = ((t124 * t143 + t146 * t149) * t141 - t143 * t158) * t136;
	t117 = t128 * t124 + t127 * t149;
	t116 = t128 * t123 - t127 * t130;
	t114 = t139 * t121 + (t123 * t146 + t130 * t143) * t136;
	t113 = t115 * t145 + t125 * t141;
	t107 = -t121 * t157 + t123 * t150 + t126 * t130;
	t105 = 0.1e1 / t112;
	t104 = t128 * t135 * t141 + t117 * t145;
	t103 = (t126 * t124 + t149 * t150) * t141 - t126 * t158;
	t100 = t108 * t145 + t141 * t147;
	t96 = 0.1e1 / (t106 * t98 ^ 2 + 0.1e1);
	t90 = 0.1e1 / t93;
	t89 = 0.1e1 / t153;
	t86 = 0.1e1 / t88;
	t85 = 0.1e1 / (t101 ^ 2 * t87 + 0.1e1);
	t84 = (-t103 * t105 + t118 * t162) * t96;
	t83 = (-t105 * t107 + t114 * t162) * t96 * t141;
	t82 = (-t100 * t105 + t113 * t162) * t96;
	t1 = [0, t84, t83, 0, t82, 0; 0, ((t117 * t141 - t128 * t158) * t86 - (-t94 * t103 + t95 * t118 + t152 * t84) * t163) * t85, (t110 * t141 * t86 - (t152 * t83 + (-t107 * t94 + t114 * t95) * t141) * t163) * t85, 0, (t102 * t86 - (-t94 * t100 + t95 * t113 + t152 * t82) * t163) * t85, 0; 0, ((t104 * t140 - t116 * t144) * t90 - (t104 * t144 + t116 * t140) * t164) * t89, ((t140 * t160 - t144 * t148) * t90 - (t140 * t148 + t144 * t160) * t164) * t89, 0, (-t140 * t90 + t144 * t164) * t89 * t101, t153 * t89;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end