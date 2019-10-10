% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR1
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
%   Wie in S6PPRRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:18
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRRR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (90->0), div. (5->0), fcn. (119->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:21
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
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:21
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (429->30), mult. (1141->67), div. (40->9), fcn. (1556->15), ass. (0->44)
	t82 = sin(pkin(13));
	t83 = sin(pkin(12));
	t86 = cos(pkin(13));
	t87 = cos(pkin(12));
	t88 = cos(pkin(7));
	t89 = cos(pkin(6));
	t97 = t87 * t89;
	t84 = sin(pkin(7));
	t85 = sin(pkin(6));
	t98 = t85 * t84;
	t101 = (-t83 * t82 + t86 * t97) * t88 - t87 * t98;
	t100 = t83 * t89;
	t99 = t84 * t89;
	t91 = cos(qJ(3));
	t96 = t88 * t91;
	t95 = t83 * t98;
	t76 = -t86 * t100 - t87 * t82;
	t77 = -t82 * t100 + t87 * t86;
	t90 = sin(qJ(3));
	t68 = t77 * t91 + (t76 * t88 + t95) * t90;
	t72 = t83 * t85 * t88 - t76 * t84;
	t81 = qJ(4) + qJ(5);
	t79 = sin(t81);
	t80 = cos(t81);
	t59 = t68 * t80 + t72 * t79;
	t57 = 0.1e1 / t59 ^ 2;
	t58 = t68 * t79 - t72 * t80;
	t93 = t58 ^ 2 * t57 + 0.1e1;
	t75 = t82 * t97 + t83 * t86;
	t71 = t90 * t99 + (t86 * t88 * t90 + t82 * t91) * t85;
	t70 = -t91 * t99 + (t82 * t90 - t86 * t96) * t85;
	t69 = 0.1e1 / t70 ^ 2;
	t67 = -t76 * t96 + t77 * t90 - t91 * t95;
	t66 = t101 * t90 + t75 * t91;
	t64 = -t101 * t91 + t75 * t90;
	t63 = atan2(-t64, t70);
	t61 = cos(t63);
	t60 = sin(t63);
	t56 = 0.1e1 / t93;
	t55 = -t60 * t64 + t61 * t70;
	t54 = 0.1e1 / t55 ^ 2;
	t52 = (-t66 / t70 + t71 * t64 * t69) / (t64 ^ 2 * t69 + 0.1e1);
	t51 = t93 * t56;
	t1 = [0, 0, t52, 0, 0, 0; 0, 0, (t68 / t55 - (-t60 * t66 + t61 * t71 + (-t60 * t70 - t61 * t64) * t52) * t67 * t54) / (t67 ^ 2 * t54 + 0.1e1), 0, 0, 0; 0, 0, (-t79 / t59 + t80 * t58 * t57) * t67 * t56, t51, t51, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:21
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (1997->44), mult. (4405->105), div. (95->9), fcn. (6043->17), ass. (0->69)
	t120 = sin(pkin(12));
	t125 = cos(pkin(7));
	t119 = sin(pkin(13));
	t123 = cos(pkin(13));
	t124 = cos(pkin(12));
	t126 = cos(pkin(6));
	t145 = t120 * t126;
	t135 = t124 * t119 + t123 * t145;
	t121 = sin(pkin(7));
	t122 = sin(pkin(6));
	t143 = t122 * t121;
	t151 = -t120 * t143 + t135 * t125;
	t128 = sin(qJ(3));
	t130 = cos(qJ(3));
	t141 = t123 * t125;
	t144 = t121 * t126;
	t110 = t128 * t144 + (t119 * t130 + t128 * t141) * t122;
	t112 = -t123 * t143 + t126 * t125;
	t118 = qJ(4) + qJ(5);
	t116 = sin(t118);
	t117 = cos(t118);
	t101 = t110 * t116 - t112 * t117;
	t140 = t124 * t126;
	t113 = t119 * t140 + t120 * t123;
	t136 = -t120 * t119 + t123 * t140;
	t132 = -t124 * t143 + t136 * t125;
	t104 = t113 * t130 + t132 * t128;
	t142 = t122 * t125;
	t133 = -t136 * t121 - t124 * t142;
	t94 = t104 * t116 - t133 * t117;
	t93 = atan2(-t94, t101);
	t90 = sin(t93);
	t91 = cos(t93);
	t84 = t91 * t101 - t90 * t94;
	t83 = 0.1e1 / t84 ^ 2;
	t114 = -t119 * t145 + t124 * t123;
	t106 = t114 * t130 - t151 * t128;
	t131 = t120 * t142 + t135 * t121;
	t97 = t106 * t116 - t131 * t117;
	t150 = t83 * t97;
	t129 = cos(qJ(6));
	t105 = t114 * t128 + t151 * t130;
	t127 = sin(qJ(6));
	t147 = t105 * t127;
	t98 = t106 * t117 + t131 * t116;
	t89 = t98 * t129 + t147;
	t87 = 0.1e1 / t89 ^ 2;
	t146 = t105 * t129;
	t88 = t98 * t127 - t146;
	t149 = t87 * t88;
	t100 = 0.1e1 / t101 ^ 2;
	t148 = t100 * t94;
	t139 = t88 ^ 2 * t87 + 0.1e1;
	t137 = -t101 * t90 - t91 * t94;
	t109 = t130 * t144 + (-t119 * t128 + t130 * t141) * t122;
	t103 = -t113 * t128 + t132 * t130;
	t102 = t110 * t117 + t112 * t116;
	t99 = 0.1e1 / t101;
	t96 = t104 * t117 + t133 * t116;
	t92 = 0.1e1 / (t94 ^ 2 * t100 + 0.1e1);
	t86 = 0.1e1 / t89;
	t85 = 0.1e1 / t139;
	t82 = 0.1e1 / t84;
	t81 = 0.1e1 / (t97 ^ 2 * t83 + 0.1e1);
	t80 = (-t103 * t99 + t109 * t148) * t92 * t116;
	t79 = (t102 * t148 - t96 * t99) * t92;
	t78 = (-t127 * t86 + t129 * t149) * t97 * t85;
	t77 = (t98 * t82 - (t91 * t102 + t137 * t79 - t90 * t96) * t150) * t81;
	t1 = [0, 0, t80, t79, t79, 0; 0, 0, (-t105 * t116 * t82 - (t137 * t80 + (-t103 * t90 + t109 * t91) * t116) * t150) * t81, t77, t77, 0; 0, 0, ((-t106 * t129 - t117 * t147) * t86 - (t106 * t127 - t117 * t146) * t149) * t85, t78, t78, t139 * t85;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end