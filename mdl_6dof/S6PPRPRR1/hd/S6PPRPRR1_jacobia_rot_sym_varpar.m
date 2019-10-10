% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPRPRR1
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
%   Wie in S6PPRPRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:08
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRPRR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRPRR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (90->0), div. (5->0), fcn. (119->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (66->0), mult. (186->0), div. (6->0), fcn. (252->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (582->31), mult. (1659->71), div. (35->9), fcn. (2273->17), ass. (0->48)
	t83 = sin(pkin(11));
	t85 = sin(pkin(6));
	t100 = t83 * t85;
	t90 = cos(pkin(6));
	t99 = t83 * t90;
	t88 = cos(pkin(11));
	t98 = t85 * t88;
	t97 = t88 * t90;
	t84 = sin(pkin(7));
	t81 = sin(pkin(13));
	t86 = cos(pkin(13));
	t92 = sin(qJ(3));
	t94 = cos(qJ(3));
	t95 = t94 * t81 + t92 * t86;
	t72 = t95 * t84;
	t89 = cos(pkin(7));
	t74 = t95 * t89;
	t82 = sin(pkin(12));
	t87 = cos(pkin(12));
	t77 = -t88 * t82 - t87 * t99;
	t78 = -t82 * t99 + t88 * t87;
	t79 = t92 * t81 - t94 * t86;
	t65 = t72 * t100 + t77 * t74 - t78 * t79;
	t69 = t89 * t100 - t77 * t84;
	t91 = sin(qJ(5));
	t93 = cos(qJ(5));
	t60 = t65 * t93 + t69 * t91;
	t58 = 0.1e1 / t60 ^ 2;
	t59 = t65 * t91 - t69 * t93;
	t96 = t59 ^ 2 * t58 + 0.1e1;
	t76 = t82 * t97 + t83 * t87;
	t75 = -t83 * t82 + t87 * t97;
	t73 = t79 * t89;
	t71 = t79 * t84;
	t68 = t90 * t71 + (t73 * t87 + t82 * t95) * t85;
	t67 = t90 * t72 + (t74 * t87 - t79 * t82) * t85;
	t66 = 0.1e1 / t68 ^ 2;
	t63 = -t71 * t100 - t77 * t73 - t78 * t95;
	t62 = t71 * t98 - t75 * t73 - t76 * t95;
	t61 = t72 * t98 - t75 * t74 + t76 * t79;
	t57 = atan2(t62, t68);
	t55 = cos(t57);
	t54 = sin(t57);
	t53 = 0.1e1 / t96;
	t52 = t54 * t62 + t55 * t68;
	t51 = 0.1e1 / t52 ^ 2;
	t49 = (t61 / t68 - t67 * t62 * t66) / (t62 ^ 2 * t66 + 0.1e1);
	t1 = [0, 0, t49, 0, 0, 0; 0, 0, (t65 / t52 + (t54 * t61 + t55 * t67 + (-t54 * t68 + t55 * t62) * t49) * t63 * t51) / (t63 ^ 2 * t51 + 0.1e1), 0, 0, 0; 0, 0, (t91 / t60 - t93 * t59 * t58) * t63 * t53, 0, t96 * t53, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:44
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (1594->47), mult. (4461->113), div. (65->9), fcn. (6121->19), ass. (0->71)
	t116 = sin(qJ(5));
	t119 = cos(qJ(5));
	t108 = sin(pkin(7));
	t112 = cos(pkin(11));
	t106 = sin(pkin(12));
	t107 = sin(pkin(11));
	t111 = cos(pkin(12));
	t114 = cos(pkin(6));
	t128 = t112 * t114;
	t123 = -t106 * t107 + t111 * t128;
	t109 = sin(pkin(6));
	t113 = cos(pkin(7));
	t129 = t109 * t113;
	t121 = -t108 * t123 - t112 * t129;
	t105 = sin(pkin(13));
	t110 = cos(pkin(13));
	t117 = sin(qJ(3));
	t120 = cos(qJ(3));
	t102 = t105 * t117 - t120 * t110;
	t130 = t109 * t112;
	t125 = t105 * t120 + t110 * t117;
	t95 = t125 * t108;
	t97 = t125 * t113;
	t99 = t106 * t128 + t107 * t111;
	t82 = -t99 * t102 + t123 * t97 - t130 * t95;
	t76 = t82 * t116 - t119 * t121;
	t91 = t114 * t95 + (-t102 * t106 + t111 * t97) * t109;
	t98 = -t108 * t109 * t111 + t113 * t114;
	t88 = t116 * t91 - t119 * t98;
	t75 = atan2(-t76, t88);
	t72 = sin(t75);
	t73 = cos(t75);
	t66 = -t72 * t76 + t73 * t88;
	t65 = 0.1e1 / t66 ^ 2;
	t131 = t107 * t114;
	t100 = -t106 * t112 - t111 * t131;
	t101 = -t106 * t131 + t111 * t112;
	t132 = t107 * t109;
	t122 = t100 * t97 - t101 * t102 + t132 * t95;
	t124 = -t100 * t108 + t107 * t129;
	t79 = t116 * t122 - t119 * t124;
	t137 = t65 * t79;
	t118 = cos(qJ(6));
	t115 = sin(qJ(6));
	t94 = t102 * t108;
	t96 = t102 * t113;
	t84 = -t100 * t96 - t101 * t125 - t132 * t94;
	t134 = t84 * t115;
	t80 = t116 * t124 + t119 * t122;
	t71 = t118 * t80 - t134;
	t69 = 0.1e1 / t71 ^ 2;
	t133 = t84 * t118;
	t70 = t80 * t115 + t133;
	t136 = t69 * t70;
	t87 = 0.1e1 / t88 ^ 2;
	t135 = t76 * t87;
	t127 = t69 * t70 ^ 2 + 0.1e1;
	t126 = -t72 * t88 - t73 * t76;
	t90 = -t114 * t94 + (-t106 * t125 - t111 * t96) * t109;
	t89 = t116 * t98 + t119 * t91;
	t86 = 0.1e1 / t88;
	t81 = -t123 * t96 - t125 * t99 + t130 * t94;
	t78 = t116 * t121 + t82 * t119;
	t74 = 0.1e1 / (t76 ^ 2 * t87 + 0.1e1);
	t68 = 0.1e1 / t71;
	t67 = 0.1e1 / t127;
	t64 = 0.1e1 / t66;
	t63 = 0.1e1 / (t65 * t79 ^ 2 + 0.1e1);
	t62 = (t135 * t90 - t81 * t86) * t74 * t116;
	t61 = (t135 * t89 - t78 * t86) * t74;
	t1 = [0, 0, t62, 0, t61, 0; 0, 0, (t84 * t116 * t64 - (t126 * t62 + (-t72 * t81 + t73 * t90) * t116) * t137) * t63, 0, (t80 * t64 - (t126 * t61 - t72 * t78 + t73 * t89) * t137) * t63, 0; 0, 0, ((-t118 * t122 + t119 * t134) * t68 - (t115 * t122 + t119 * t133) * t136) * t67, 0, (-t115 * t68 + t118 * t136) * t79 * t67, t127 * t67;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end