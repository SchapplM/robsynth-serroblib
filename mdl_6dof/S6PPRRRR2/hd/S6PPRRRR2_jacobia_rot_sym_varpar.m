% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR2
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
%   Wie in S6PPRRRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:20
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRRR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (90->0), div. (5->0), fcn. (119->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (333->29), mult. (995->67), div. (35->9), fcn. (1360->15), ass. (0->42)
	t64 = sin(pkin(13));
	t65 = sin(pkin(12));
	t68 = cos(pkin(13));
	t69 = cos(pkin(12));
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
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (1048->42), mult. (3005->102), div. (65->9), fcn. (4119->17), ass. (0->66)
	t101 = cos(pkin(7));
	t100 = cos(pkin(12));
	t102 = cos(pkin(6));
	t96 = sin(pkin(12));
	t123 = t102 * t96;
	t95 = sin(pkin(13));
	t99 = cos(pkin(13));
	t115 = t100 * t95 + t99 * t123;
	t97 = sin(pkin(7));
	t98 = sin(pkin(6));
	t124 = t98 * t97;
	t129 = t115 * t101 - t96 * t124;
	t104 = sin(qJ(4));
	t107 = cos(qJ(4));
	t119 = t100 * t102;
	t114 = t99 * t119 - t96 * t95;
	t120 = t98 * t101;
	t109 = -t100 * t120 - t114 * t97;
	t105 = sin(qJ(3));
	t108 = cos(qJ(3));
	t110 = -t100 * t124 + t114 * t101;
	t92 = t95 * t119 + t96 * t99;
	t79 = t110 * t105 + t92 * t108;
	t73 = t79 * t104 - t109 * t107;
	t113 = t102 * t97 + t99 * t120;
	t125 = t95 * t98;
	t89 = t113 * t105 + t108 * t125;
	t91 = t102 * t101 - t99 * t124;
	t84 = t89 * t104 - t91 * t107;
	t72 = atan2(-t73, t84);
	t69 = sin(t72);
	t70 = cos(t72);
	t63 = -t69 * t73 + t70 * t84;
	t62 = 0.1e1 / t63 ^ 2;
	t111 = t115 * t97 + t96 * t120;
	t93 = t100 * t99 - t95 * t123;
	t81 = -t129 * t105 + t93 * t108;
	t76 = t81 * t104 - t111 * t107;
	t128 = t62 * t76;
	t106 = cos(qJ(5));
	t103 = sin(qJ(5));
	t80 = t93 * t105 + t129 * t108;
	t122 = t80 * t103;
	t77 = t111 * t104 + t81 * t107;
	t68 = t77 * t106 + t122;
	t66 = 0.1e1 / t68 ^ 2;
	t121 = t80 * t106;
	t67 = t77 * t103 - t121;
	t127 = t66 * t67;
	t83 = 0.1e1 / t84 ^ 2;
	t126 = t73 * t83;
	t117 = t67 ^ 2 * t66 + 0.1e1;
	t116 = -t69 * t84 - t70 * t73;
	t88 = -t105 * t125 + t113 * t108;
	t85 = t91 * t104 + t89 * t107;
	t82 = 0.1e1 / t84;
	t78 = -t92 * t105 + t110 * t108;
	t75 = t109 * t104 + t79 * t107;
	t71 = 0.1e1 / (t73 ^ 2 * t83 + 0.1e1);
	t65 = 0.1e1 / t68;
	t64 = 0.1e1 / t117;
	t61 = 0.1e1 / t63;
	t60 = 0.1e1 / (t76 ^ 2 * t62 + 0.1e1);
	t59 = (t88 * t126 - t78 * t82) * t71 * t104;
	t58 = (t85 * t126 - t75 * t82) * t71;
	t1 = [0, 0, t59, t58, 0, 0; 0, 0, (-t80 * t104 * t61 - (t116 * t59 + (-t69 * t78 + t70 * t88) * t104) * t128) * t60, (t77 * t61 - (t116 * t58 - t69 * t75 + t70 * t85) * t128) * t60, 0, 0; 0, 0, ((-t81 * t106 - t107 * t122) * t65 - (t81 * t103 - t107 * t121) * t127) * t64, (-t103 * t65 + t106 * t127) * t76 * t64, t117 * t64, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (1197->44), mult. (3249->105), div. (70->9), fcn. (4448->17), ass. (0->68)
	t114 = sin(pkin(12));
	t119 = cos(pkin(7));
	t113 = sin(pkin(13));
	t117 = cos(pkin(13));
	t118 = cos(pkin(12));
	t120 = cos(pkin(6));
	t139 = t114 * t120;
	t129 = t118 * t113 + t117 * t139;
	t115 = sin(pkin(7));
	t116 = sin(pkin(6));
	t137 = t116 * t115;
	t145 = -t114 * t137 + t129 * t119;
	t121 = sin(qJ(4));
	t123 = cos(qJ(4));
	t134 = t118 * t120;
	t130 = -t114 * t113 + t117 * t134;
	t136 = t116 * t119;
	t127 = -t130 * t115 - t118 * t136;
	t107 = t113 * t134 + t114 * t117;
	t122 = sin(qJ(3));
	t124 = cos(qJ(3));
	t126 = -t118 * t137 + t130 * t119;
	t94 = t107 * t124 + t126 * t122;
	t88 = t94 * t121 - t127 * t123;
	t135 = t117 * t119;
	t138 = t115 * t120;
	t104 = t122 * t138 + (t113 * t124 + t122 * t135) * t116;
	t106 = -t117 * t137 + t120 * t119;
	t99 = t104 * t121 - t106 * t123;
	t87 = atan2(-t88, t99);
	t84 = sin(t87);
	t85 = cos(t87);
	t78 = -t84 * t88 + t85 * t99;
	t77 = 0.1e1 / t78 ^ 2;
	t125 = t114 * t136 + t129 * t115;
	t108 = -t113 * t139 + t118 * t117;
	t96 = t108 * t124 - t145 * t122;
	t91 = t96 * t121 - t125 * t123;
	t144 = t77 * t91;
	t112 = qJ(5) + qJ(6);
	t111 = cos(t112);
	t110 = sin(t112);
	t95 = t108 * t122 + t145 * t124;
	t141 = t95 * t110;
	t92 = t125 * t121 + t96 * t123;
	t83 = t92 * t111 + t141;
	t81 = 0.1e1 / t83 ^ 2;
	t140 = t95 * t111;
	t82 = t92 * t110 - t140;
	t143 = t81 * t82;
	t98 = 0.1e1 / t99 ^ 2;
	t142 = t88 * t98;
	t133 = t82 ^ 2 * t81 + 0.1e1;
	t131 = -t84 * t99 - t85 * t88;
	t103 = t124 * t138 + (-t113 * t122 + t124 * t135) * t116;
	t100 = t104 * t123 + t106 * t121;
	t97 = 0.1e1 / t99;
	t93 = -t107 * t122 + t126 * t124;
	t90 = t127 * t121 + t94 * t123;
	t86 = 0.1e1 / (t88 ^ 2 * t98 + 0.1e1);
	t80 = 0.1e1 / t83;
	t79 = 0.1e1 / t133;
	t76 = 0.1e1 / t78;
	t75 = 0.1e1 / (t91 ^ 2 * t77 + 0.1e1);
	t74 = (t103 * t142 - t93 * t97) * t86 * t121;
	t73 = (t100 * t142 - t90 * t97) * t86;
	t72 = t133 * t79;
	t1 = [0, 0, t74, t73, 0, 0; 0, 0, (-t95 * t121 * t76 - (t131 * t74 + (t103 * t85 - t84 * t93) * t121) * t144) * t75, (t92 * t76 - (t85 * t100 + t131 * t73 - t84 * t90) * t144) * t75, 0, 0; 0, 0, ((-t96 * t111 - t123 * t141) * t80 - (t96 * t110 - t123 * t140) * t143) * t79, (-t110 * t80 + t111 * t143) * t91 * t79, t72, t72;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end