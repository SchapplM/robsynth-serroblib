% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRPR2
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
%   Wie in S6PPRRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRPR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (90->0), div. (5->0), fcn. (119->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (333->29), mult. (995->67), div. (35->9), fcn. (1360->15), ass. (0->42)
	t64 = sin(pkin(12));
	t65 = sin(pkin(11));
	t68 = cos(pkin(12));
	t69 = cos(pkin(11));
	t70 = cos(pkin(7));
	t71 = cos(pkin(6));
	t81 = t69 * t71;
	t66 = sin(pkin(7));
	t67 = sin(pkin(6));
	t82 = t67 * t66;
	t85 = (-t65 * t64 + t68 * t81) * t70 - t69 * t82;
	t84 = t65 * t71;
	t83 = t66 * t71;
	t75 = cos(qJ(3));
	t80 = t70 * t75;
	t79 = t65 * t82;
	t61 = -t69 * t64 - t68 * t84;
	t62 = -t64 * t84 + t69 * t68;
	t73 = sin(qJ(3));
	t53 = t62 * t75 + (t61 * t70 + t79) * t73;
	t57 = t65 * t67 * t70 - t61 * t66;
	t72 = sin(qJ(4));
	t74 = cos(qJ(4));
	t44 = t53 * t74 + t57 * t72;
	t42 = 0.1e1 / t44 ^ 2;
	t43 = t53 * t72 - t57 * t74;
	t77 = t43 ^ 2 * t42 + 0.1e1;
	t60 = t64 * t81 + t65 * t68;
	t56 = t73 * t83 + (t68 * t70 * t73 + t64 * t75) * t67;
	t55 = -t75 * t83 + (t64 * t73 - t68 * t80) * t67;
	t54 = 0.1e1 / t55 ^ 2;
	t52 = -t61 * t80 + t62 * t73 - t75 * t79;
	t51 = t60 * t75 + t85 * t73;
	t49 = t60 * t73 - t85 * t75;
	t48 = atan2(-t49, t55);
	t46 = cos(t48);
	t45 = sin(t48);
	t41 = 0.1e1 / t77;
	t40 = -t45 * t49 + t46 * t55;
	t39 = 0.1e1 / t40 ^ 2;
	t37 = (-t51 / t55 + t56 * t49 * t54) / (t49 ^ 2 * t54 + 0.1e1);
	t1 = [0, 0, t37, 0, 0, 0; 0, 0, (t53 / t40 - (-t45 * t51 + t46 * t56 + (-t45 * t55 - t46 * t49) * t37) * t52 * t39) / (t52 ^ 2 * t39 + 0.1e1), 0, 0, 0; 0, 0, (-t72 / t44 + t74 * t43 * t42) * t52 * t41, t77 * t41, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (826->40), mult. (2382->97), div. (57->9), fcn. (3277->15), ass. (0->58)
	t81 = sin(pkin(12));
	t82 = sin(pkin(11));
	t85 = cos(pkin(12));
	t86 = cos(pkin(11));
	t88 = cos(pkin(6));
	t99 = t86 * t88;
	t77 = t81 * t99 + t82 * t85;
	t90 = sin(qJ(3));
	t92 = cos(qJ(3));
	t83 = sin(pkin(7));
	t84 = sin(pkin(6));
	t101 = t84 * t83;
	t87 = cos(pkin(7));
	t95 = -t82 * t81 + t85 * t99;
	t93 = -t86 * t101 + t95 * t87;
	t64 = t77 * t92 + t93 * t90;
	t89 = sin(qJ(4));
	t91 = cos(qJ(4));
	t100 = t84 * t87;
	t94 = -t86 * t100 - t95 * t83;
	t56 = t64 * t89 - t94 * t91;
	t102 = t83 * t88;
	t73 = t90 * t102 + (t85 * t87 * t90 + t81 * t92) * t84;
	t76 = -t85 * t101 + t88 * t87;
	t69 = t73 * t89 - t76 * t91;
	t55 = atan2(-t56, t69);
	t52 = sin(t55);
	t53 = cos(t55);
	t50 = -t52 * t56 + t53 * t69;
	t49 = 0.1e1 / t50 ^ 2;
	t103 = t82 * t88;
	t78 = -t85 * t103 - t86 * t81;
	t79 = -t81 * t103 + t86 * t85;
	t97 = t82 * t101;
	t66 = t79 * t92 + (t78 * t87 + t97) * t90;
	t74 = t82 * t100 - t78 * t83;
	t59 = t66 * t89 - t74 * t91;
	t105 = t49 * t59;
	t68 = 0.1e1 / t69 ^ 2;
	t104 = t56 * t68;
	t98 = t87 * t92;
	t96 = -t52 * t69 - t53 * t56;
	t72 = t92 * t102 + (-t81 * t90 + t85 * t98) * t84;
	t70 = t73 * t91 + t76 * t89;
	t67 = 0.1e1 / t69;
	t65 = -t78 * t98 + t79 * t90 - t92 * t97;
	t63 = -t77 * t90 + t93 * t92;
	t62 = 0.1e1 / t65 ^ 2;
	t61 = 0.1e1 / t65;
	t60 = t66 * t91 + t74 * t89;
	t58 = t64 * t91 + t94 * t89;
	t54 = 0.1e1 / (t56 ^ 2 * t68 + 0.1e1);
	t51 = 0.1e1 / (t60 ^ 2 * t62 + 0.1e1);
	t48 = 0.1e1 / t50;
	t47 = 0.1e1 / (t59 ^ 2 * t49 + 0.1e1);
	t46 = (t72 * t104 - t63 * t67) * t89 * t54;
	t45 = (t70 * t104 - t58 * t67) * t54;
	t1 = [0, 0, t46, t45, 0, 0; 0, 0, (-t65 * t89 * t48 - ((-t52 * t63 + t53 * t72) * t89 + t96 * t46) * t105) * t47, (t60 * t48 - (t96 * t45 - t52 * t58 + t53 * t70) * t105) * t47, 0, 0; 0, 0, (-t60 * t62 * t66 - t61 * t65 * t91) * t51, -t59 * t61 * t51, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (1048->43), mult. (3005->105), div. (65->9), fcn. (4119->17), ass. (0->66)
	t102 = sin(pkin(11));
	t107 = cos(pkin(7));
	t101 = sin(pkin(12));
	t105 = cos(pkin(12));
	t106 = cos(pkin(11));
	t108 = cos(pkin(6));
	t126 = t102 * t108;
	t116 = t106 * t101 + t105 * t126;
	t103 = sin(pkin(7));
	t104 = sin(pkin(6));
	t124 = t104 * t103;
	t132 = -t102 * t124 + t116 * t107;
	t110 = sin(qJ(4));
	t113 = cos(qJ(4));
	t111 = sin(qJ(3));
	t114 = cos(qJ(3));
	t121 = t106 * t108;
	t97 = -t102 * t101 + t105 * t121;
	t117 = -t106 * t124 + t107 * t97;
	t98 = t101 * t121 + t102 * t105;
	t84 = t117 * t111 + t98 * t114;
	t123 = t104 * t107;
	t93 = -t97 * t103 - t106 * t123;
	t79 = t93 * t110 + t84 * t113;
	t122 = t105 * t107;
	t125 = t103 * t108;
	t92 = t111 * t125 + (t101 * t114 + t111 * t122) * t104;
	t96 = -t105 * t124 + t108 * t107;
	t90 = t96 * t110 + t92 * t113;
	t77 = atan2(-t79, t90);
	t74 = sin(t77);
	t75 = cos(t77);
	t68 = -t74 * t79 + t75 * t90;
	t67 = 0.1e1 / t68 ^ 2;
	t99 = -t101 * t126 + t106 * t105;
	t86 = -t132 * t111 + t99 * t114;
	t94 = t102 * t123 + t116 * t103;
	t82 = t94 * t110 + t86 * t113;
	t131 = t67 * t82;
	t109 = sin(qJ(6));
	t112 = cos(qJ(6));
	t85 = t99 * t111 + t132 * t114;
	t127 = t85 * t112;
	t81 = t86 * t110 - t94 * t113;
	t73 = t81 * t109 + t127;
	t71 = 0.1e1 / t73 ^ 2;
	t128 = t85 * t109;
	t72 = -t81 * t112 + t128;
	t130 = t71 * t72;
	t88 = 0.1e1 / t90 ^ 2;
	t129 = t79 * t88;
	t120 = t72 ^ 2 * t71 + 0.1e1;
	t118 = -t74 * t90 - t75 * t79;
	t91 = t114 * t125 + (-t101 * t111 + t114 * t122) * t104;
	t89 = -t92 * t110 + t96 * t113;
	t87 = 0.1e1 / t90;
	t83 = -t98 * t111 + t117 * t114;
	t78 = t84 * t110 - t93 * t113;
	t76 = 0.1e1 / (t79 ^ 2 * t88 + 0.1e1);
	t70 = 0.1e1 / t73;
	t69 = 0.1e1 / t120;
	t66 = 0.1e1 / t68;
	t65 = 0.1e1 / (t82 ^ 2 * t67 + 0.1e1);
	t64 = (t91 * t129 - t83 * t87) * t76 * t113;
	t63 = (t89 * t129 + t78 * t87) * t76;
	t1 = [0, 0, t64, t63, 0, 0; 0, 0, (-t85 * t113 * t66 - (t118 * t64 + (-t74 * t83 + t75 * t91) * t113) * t131) * t65, (-t81 * t66 - (t118 * t63 + t74 * t78 + t75 * t89) * t131) * t65, 0, 0; 0, 0, ((t86 * t109 + t110 * t127) * t70 - (-t110 * t128 + t86 * t112) * t130) * t69, (-t109 * t130 - t112 * t70) * t82 * t69, 0, t120 * t69;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end