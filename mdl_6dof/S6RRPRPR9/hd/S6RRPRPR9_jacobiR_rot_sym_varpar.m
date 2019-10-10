% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:19
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR9_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR9_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t50 = sin(qJ(2));
	t51 = sin(qJ(1));
	t57 = t51 * t50;
	t52 = cos(qJ(2));
	t56 = t51 * t52;
	t53 = cos(qJ(1));
	t55 = t53 * t50;
	t54 = t53 * t52;
	t49 = cos(pkin(6));
	t48 = sin(pkin(6));
	t47 = -t49 * t57 + t54;
	t46 = -t49 * t56 - t55;
	t45 = -t49 * t55 - t56;
	t44 = -t49 * t54 + t57;
	t1 = [t45, t46, 0, 0, 0, 0; t47, -t44, 0, 0, 0, 0; 0, t48 * t52, 0, 0, 0, 0; t44, -t47, 0, 0, 0, 0; t46, t45, 0, 0, 0, 0; 0, -t48 * t50, 0, 0, 0, 0; t53 * t48, 0, 0, 0, 0, 0; t51 * t48, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (19->11), mult. (57->26), div. (0->0), fcn. (88->8), ass. (0->20)
	t66 = sin(pkin(6));
	t70 = sin(qJ(1));
	t79 = t66 * t70;
	t71 = cos(qJ(2));
	t78 = t66 * t71;
	t72 = cos(qJ(1));
	t77 = t66 * t72;
	t69 = sin(qJ(2));
	t76 = t70 * t69;
	t75 = t70 * t71;
	t74 = t72 * t69;
	t73 = t72 * t71;
	t68 = cos(pkin(6));
	t67 = cos(pkin(11));
	t65 = sin(pkin(11));
	t64 = -t68 * t76 + t73;
	t63 = t68 * t75 + t74;
	t62 = t68 * t74 + t75;
	t61 = t68 * t73 - t76;
	t1 = [-t62 * t67 + t65 * t77, -t63 * t67, 0, 0, 0, 0; t64 * t67 + t65 * t79, t61 * t67, 0, 0, 0, 0; 0, t67 * t78, 0, 0, 0, 0; t62 * t65 + t67 * t77, t63 * t65, 0, 0, 0, 0; -t64 * t65 + t67 * t79, -t61 * t65, 0, 0, 0, 0; 0, -t65 * t78, 0, 0, 0, 0; t61, t64, 0, 0, 0, 0; t63, t62, 0, 0, 0, 0; 0, t66 * t69, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (55->16), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->26)
	t90 = sin(pkin(6));
	t92 = sin(qJ(2));
	t105 = t90 * t92;
	t93 = sin(qJ(1));
	t104 = t90 * t93;
	t94 = cos(qJ(2));
	t103 = t90 * t94;
	t95 = cos(qJ(1));
	t102 = t90 * t95;
	t101 = t93 * t92;
	t100 = t93 * t94;
	t99 = t95 * t92;
	t98 = t95 * t94;
	t91 = cos(pkin(6));
	t83 = t91 * t99 + t100;
	t89 = pkin(11) + qJ(4);
	t87 = sin(t89);
	t88 = cos(t89);
	t97 = t87 * t102 - t83 * t88;
	t96 = t88 * t102 + t83 * t87;
	t85 = -t91 * t101 + t98;
	t84 = t91 * t100 + t99;
	t82 = t91 * t98 - t101;
	t81 = t87 * t104 + t85 * t88;
	t80 = t88 * t104 - t85 * t87;
	t1 = [t97, -t84 * t88, 0, t80, 0, 0; t81, t82 * t88, 0, -t96, 0, 0; 0, t88 * t103, 0, -t87 * t105 + t91 * t88, 0, 0; t96, t84 * t87, 0, -t81, 0, 0; t80, -t82 * t87, 0, t97, 0, 0; 0, -t87 * t103, 0, -t88 * t105 - t91 * t87, 0, 0; t82, t85, 0, 0, 0, 0; t84, t83, 0, 0, 0, 0; 0, t105, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:18:59
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (93->26), mult. (163->58), div. (0->0), fcn. (238->10), ass. (0->31)
	t120 = pkin(11) + qJ(4);
	t119 = cos(t120);
	t121 = sin(pkin(12));
	t138 = t119 * t121;
	t123 = cos(pkin(12));
	t137 = t119 * t123;
	t127 = cos(qJ(2));
	t136 = t119 * t127;
	t122 = sin(pkin(6));
	t125 = sin(qJ(2));
	t135 = t122 * t125;
	t126 = sin(qJ(1));
	t134 = t122 * t126;
	t128 = cos(qJ(1));
	t133 = t122 * t128;
	t132 = t126 * t125;
	t131 = t126 * t127;
	t130 = t128 * t125;
	t129 = t128 * t127;
	t124 = cos(pkin(6));
	t114 = t124 * t130 + t131;
	t118 = sin(t120);
	t108 = -t114 * t118 - t119 * t133;
	t109 = -t114 * t119 + t118 * t133;
	t116 = -t124 * t132 + t129;
	t115 = t124 * t131 + t130;
	t113 = t124 * t129 - t132;
	t112 = -t118 * t135 + t124 * t119;
	t111 = t116 * t119 + t118 * t134;
	t110 = t116 * t118 - t119 * t134;
	t1 = [t109 * t123 + t113 * t121, -t115 * t137 + t116 * t121, 0, -t110 * t123, 0, 0; t111 * t123 + t115 * t121, t113 * t137 + t114 * t121, 0, t108 * t123, 0, 0; 0, (t121 * t125 + t123 * t136) * t122, 0, t112 * t123, 0, 0; -t109 * t121 + t113 * t123, t115 * t138 + t116 * t123, 0, t110 * t121, 0, 0; -t111 * t121 + t115 * t123, -t113 * t138 + t114 * t123, 0, -t108 * t121, 0, 0; 0, (-t121 * t136 + t123 * t125) * t122, 0, -t112 * t121, 0, 0; t108, -t115 * t118, 0, t111, 0, 0; t110, t113 * t118, 0, -t109, 0, 0; 0, t122 * t127 * t118, 0, t124 * t118 + t119 * t135, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:18:59
	% EndTime: 2019-10-10 10:19:00
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (163->32), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->38)
	t146 = cos(pkin(6));
	t147 = sin(qJ(2));
	t150 = cos(qJ(1));
	t152 = t150 * t147;
	t148 = sin(qJ(1));
	t149 = cos(qJ(2));
	t153 = t148 * t149;
	t134 = t146 * t152 + t153;
	t144 = pkin(11) + qJ(4);
	t140 = sin(t144);
	t142 = cos(t144);
	t145 = sin(pkin(6));
	t155 = t145 * t150;
	t128 = -t134 * t142 + t140 * t155;
	t151 = t150 * t149;
	t154 = t148 * t147;
	t133 = -t146 * t151 + t154;
	t143 = pkin(12) + qJ(6);
	t139 = sin(t143);
	t141 = cos(t143);
	t165 = t128 * t139 + t133 * t141;
	t164 = t128 * t141 - t133 * t139;
	t161 = t139 * t142;
	t160 = t141 * t142;
	t159 = t142 * t149;
	t158 = t145 * t147;
	t157 = t145 * t148;
	t156 = t145 * t149;
	t126 = -t134 * t140 - t142 * t155;
	t136 = -t146 * t154 + t151;
	t135 = t146 * t153 + t152;
	t132 = t146 * t140 + t142 * t158;
	t131 = -t140 * t158 + t146 * t142;
	t130 = t136 * t142 + t140 * t157;
	t129 = t136 * t140 - t142 * t157;
	t125 = t130 * t141 + t135 * t139;
	t124 = -t130 * t139 + t135 * t141;
	t1 = [t164, -t135 * t160 + t136 * t139, 0, -t129 * t141, 0, t124; t125, -t133 * t160 + t134 * t139, 0, t126 * t141, 0, t165; 0, (t139 * t147 + t141 * t159) * t145, 0, t131 * t141, 0, -t132 * t139 - t141 * t156; -t165, t135 * t161 + t136 * t141, 0, t129 * t139, 0, -t125; t124, t133 * t161 + t134 * t141, 0, -t126 * t139, 0, t164; 0, (-t139 * t159 + t141 * t147) * t145, 0, -t131 * t139, 0, -t132 * t141 + t139 * t156; t126, -t135 * t140, 0, t130, 0, 0; t129, -t133 * t140, 0, -t128, 0, 0; 0, t140 * t156, 0, t132, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end