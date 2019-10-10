% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR4
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
%   Wie in S6PRRRRR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:19
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRR4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (179->22), mult. (525->56), div. (35->9), fcn. (736->13), ass. (0->37)
	t45 = sin(pkin(13));
	t48 = cos(pkin(13));
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
	t32 = t51 * t55 + t62;
	t30 = 0.1e1 / t32 ^ 2;
	t31 = -t53 * t55 + t63;
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
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (607->40), mult. (1825->97), div. (65->9), fcn. (2498->15), ass. (0->60)
	t87 = sin(pkin(7));
	t88 = sin(pkin(6));
	t109 = t88 * t87;
	t89 = cos(pkin(13));
	t90 = cos(pkin(7));
	t91 = cos(pkin(6));
	t97 = cos(qJ(2));
	t106 = t91 * t97;
	t86 = sin(pkin(13));
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
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (721->41), mult. (1971->97), div. (70->9), fcn. (2694->15), ass. (0->62)
	t105 = cos(pkin(13));
	t106 = cos(pkin(7));
	t102 = sin(pkin(13));
	t109 = sin(qJ(2));
	t107 = cos(pkin(6));
	t111 = cos(qJ(2));
	t120 = t107 * t111;
	t113 = -t102 * t109 + t105 * t120;
	t103 = sin(pkin(7));
	t104 = sin(pkin(6));
	t123 = t104 * t103;
	t130 = -t105 * t123 + t113 * t106;
	t108 = sin(qJ(3));
	t110 = cos(qJ(3));
	t121 = t107 * t109;
	t94 = t102 * t111 + t105 * t121;
	t79 = t94 * t108 - t130 * t110;
	t122 = t106 * t110;
	t124 = t103 * t107;
	t88 = -t110 * t124 + (t108 * t109 - t111 * t122) * t104;
	t78 = atan2(-t79, t88);
	t75 = sin(t78);
	t76 = cos(t78);
	t69 = -t75 * t79 + t76 * t88;
	t68 = 0.1e1 / t69 ^ 2;
	t116 = t102 * t123;
	t96 = -t102 * t121 + t105 * t111;
	t125 = t96 * t108;
	t95 = -t102 * t120 - t105 * t109;
	t82 = -t110 * t116 - t95 * t122 + t125;
	t129 = t68 * t82;
	t101 = qJ(4) + qJ(5);
	t100 = cos(t101);
	t83 = t96 * t110 + (t106 * t95 + t116) * t108;
	t90 = t102 * t104 * t106 - t95 * t103;
	t99 = sin(t101);
	t74 = t83 * t100 + t90 * t99;
	t72 = 0.1e1 / t74 ^ 2;
	t73 = -t90 * t100 + t83 * t99;
	t128 = t72 * t73;
	t87 = 0.1e1 / t88 ^ 2;
	t127 = t79 * t87;
	t126 = t103 * t96;
	t119 = t108 * t111;
	t118 = t109 * t110;
	t117 = t73 ^ 2 * t72 + 0.1e1;
	t114 = -t75 * t88 - t76 * t79;
	t93 = (t106 * t118 + t119) * t104;
	t89 = t108 * t124 + (t106 * t119 + t118) * t104;
	t86 = 0.1e1 / t88;
	t85 = -t106 * t125 + t95 * t110;
	t84 = t113 * t108 + t94 * t122;
	t81 = t130 * t108 + t94 * t110;
	t77 = 0.1e1 / (t79 ^ 2 * t87 + 0.1e1);
	t71 = 0.1e1 / t74;
	t70 = 0.1e1 / t117;
	t67 = 0.1e1 / t69;
	t66 = 0.1e1 / (t82 ^ 2 * t68 + 0.1e1);
	t65 = (t93 * t127 - t84 * t86) * t77;
	t64 = (t89 * t127 - t81 * t86) * t77;
	t63 = t117 * t70;
	t1 = [0, t65, t64, 0, 0, 0; 0, ((t95 * t108 + t96 * t122) * t67 - (t114 * t65 - t75 * t84 + t76 * t93) * t129) * t66, (t83 * t67 - (t114 * t64 - t75 * t81 + t76 * t89) * t129) * t66, 0, 0, 0; 0, ((-t100 * t126 + t85 * t99) * t71 - (t85 * t100 + t99 * t126) * t128) * t70, (t100 * t128 - t99 * t71) * t82 * t70, t63, t63, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:24
	% DurationCPUTime: 0.60s
	% Computational Cost: add. (2625->62), mult. (5778->144), div. (125->9), fcn. (7926->17), ass. (0->78)
	t144 = sin(pkin(13));
	t147 = cos(pkin(13));
	t155 = cos(qJ(2));
	t149 = cos(pkin(6));
	t152 = sin(qJ(2));
	t165 = t149 * t152;
	t136 = t144 * t155 + t147 * t165;
	t151 = sin(qJ(3));
	t154 = cos(qJ(3));
	t148 = cos(pkin(7));
	t164 = t149 * t155;
	t158 = -t144 * t152 + t147 * t164;
	t145 = sin(pkin(7));
	t146 = sin(pkin(6));
	t169 = t146 * t145;
	t156 = -t147 * t169 + t158 * t148;
	t124 = t136 * t154 + t156 * t151;
	t143 = qJ(4) + qJ(5);
	t141 = sin(t143);
	t142 = cos(t143);
	t168 = t146 * t148;
	t157 = -t158 * t145 - t147 * t168;
	t112 = t124 * t141 - t157 * t142;
	t167 = t148 * t151;
	t170 = t145 * t149;
	t133 = t151 * t170 + (t152 * t154 + t155 * t167) * t146;
	t135 = t149 * t148 - t155 * t169;
	t121 = t133 * t141 - t135 * t142;
	t111 = atan2(-t112, t121);
	t108 = sin(t111);
	t109 = cos(t111);
	t102 = -t108 * t112 + t109 * t121;
	t101 = 0.1e1 / t102 ^ 2;
	t137 = -t144 * t164 - t147 * t152;
	t138 = -t144 * t165 + t147 * t155;
	t161 = t144 * t169;
	t126 = t138 * t154 + (t137 * t148 + t161) * t151;
	t159 = -t137 * t145 + t144 * t168;
	t115 = t126 * t141 - t159 * t142;
	t176 = t101 * t115;
	t116 = t126 * t142 + t159 * t141;
	t153 = cos(qJ(6));
	t166 = t148 * t154;
	t125 = -t137 * t166 + t138 * t151 - t154 * t161;
	t150 = sin(qJ(6));
	t173 = t125 * t150;
	t107 = t116 * t153 + t173;
	t105 = 0.1e1 / t107 ^ 2;
	t172 = t125 * t153;
	t106 = t116 * t150 - t172;
	t175 = t105 * t106;
	t120 = 0.1e1 / t121 ^ 2;
	t174 = t112 * t120;
	t171 = t145 * t142;
	t163 = t151 * t152;
	t162 = t154 * t155;
	t160 = t106 ^ 2 * t105 + 0.1e1;
	t132 = t154 * t170 + (t148 * t162 - t163) * t146;
	t129 = ((-t148 * t163 + t162) * t141 - t152 * t171) * t146;
	t128 = t137 * t154 - t138 * t167;
	t127 = t137 * t151 + t138 * t166;
	t123 = -t136 * t151 + t156 * t154;
	t122 = t133 * t142 + t135 * t141;
	t119 = 0.1e1 / t121;
	t118 = t138 * t145 * t141 + t128 * t142;
	t117 = (-t136 * t167 + t158 * t154) * t141 - t136 * t171;
	t114 = t124 * t142 + t157 * t141;
	t110 = 0.1e1 / (t112 ^ 2 * t120 + 0.1e1);
	t104 = 0.1e1 / t107;
	t103 = 0.1e1 / t160;
	t100 = 0.1e1 / t102;
	t99 = 0.1e1 / (t115 ^ 2 * t101 + 0.1e1);
	t98 = (-t119 * t123 + t132 * t174) * t141 * t110;
	t97 = (-t117 * t119 + t129 * t174) * t110;
	t96 = (-t114 * t119 + t122 * t174) * t110;
	t95 = (-t150 * t104 + t153 * t175) * t115 * t103;
	t94 = (t116 * t100 - ((-t112 * t96 + t122) * t109 + (-t121 * t96 - t114) * t108) * t176) * t99;
	t1 = [0, t97, t98, t96, t96, 0; 0, ((t128 * t141 - t138 * t171) * t100 - ((-t112 * t97 + t129) * t109 + (-t121 * t97 - t117) * t108) * t176) * t99, (-t125 * t141 * t100 - ((-t112 * t98 + t132 * t141) * t109 + (-t121 * t98 - t123 * t141) * t108) * t176) * t99, t94, t94, 0; 0, ((t118 * t150 - t127 * t153) * t104 - (t118 * t153 + t127 * t150) * t175) * t103, ((-t126 * t153 - t142 * t173) * t104 - (t126 * t150 - t142 * t172) * t175) * t103, t95, t95, t160 * t103;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end