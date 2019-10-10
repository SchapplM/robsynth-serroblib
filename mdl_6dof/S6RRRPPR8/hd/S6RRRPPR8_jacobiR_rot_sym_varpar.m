% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
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
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (29->15), mult. (87->31), div. (0->0), fcn. (134->8), ass. (0->25)
	t83 = sin(pkin(6));
	t86 = sin(qJ(2));
	t100 = t83 * t86;
	t88 = cos(qJ(3));
	t99 = t83 * t88;
	t89 = cos(qJ(2));
	t98 = t83 * t89;
	t90 = cos(qJ(1));
	t97 = t83 * t90;
	t87 = sin(qJ(1));
	t96 = t87 * t86;
	t95 = t87 * t89;
	t94 = t90 * t86;
	t93 = t90 * t89;
	t84 = cos(pkin(6));
	t79 = t84 * t94 + t95;
	t85 = sin(qJ(3));
	t92 = -t79 * t88 + t85 * t97;
	t91 = t79 * t85 + t88 * t97;
	t81 = -t84 * t96 + t93;
	t80 = t84 * t95 + t94;
	t78 = t84 * t93 - t96;
	t77 = t87 * t83 * t85 + t81 * t88;
	t76 = -t81 * t85 + t87 * t99;
	t1 = [t92, -t80 * t88, t76, 0, 0, 0; t77, t78 * t88, -t91, 0, 0, 0; 0, t88 * t98, -t85 * t100 + t84 * t88, 0, 0, 0; t91, t80 * t85, -t77, 0, 0, 0; t76, -t78 * t85, t92, 0, 0, 0; 0, -t85 * t98, -t84 * t85 - t86 * t99, 0, 0, 0; t78, t81, 0, 0, 0, 0; t80, t79, 0, 0, 0, 0; 0, t100, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (26->14), mult. (87->31), div. (0->0), fcn. (134->8), ass. (0->25)
	t102 = sin(pkin(6));
	t105 = sin(qJ(2));
	t119 = t102 * t105;
	t107 = cos(qJ(3));
	t118 = t102 * t107;
	t108 = cos(qJ(2));
	t117 = t102 * t108;
	t109 = cos(qJ(1));
	t116 = t102 * t109;
	t106 = sin(qJ(1));
	t115 = t106 * t105;
	t114 = t106 * t108;
	t113 = t109 * t105;
	t112 = t109 * t108;
	t104 = sin(qJ(3));
	t103 = cos(pkin(6));
	t99 = t103 * t113 + t114;
	t111 = -t99 * t104 - t107 * t116;
	t110 = t104 * t116 - t99 * t107;
	t101 = -t103 * t115 + t112;
	t100 = t103 * t114 + t113;
	t98 = t103 * t112 - t115;
	t97 = t106 * t102 * t104 + t101 * t107;
	t96 = t101 * t104 - t106 * t118;
	t1 = [t110, -t100 * t107, -t96, 0, 0, 0; t97, t98 * t107, t111, 0, 0, 0; 0, t107 * t117, t103 * t107 - t104 * t119, 0, 0, 0; t98, t101, 0, 0, 0, 0; t100, t99, 0, 0, 0, 0; 0, t119, 0, 0, 0, 0; t111, -t100 * t104, t97, 0, 0, 0; t96, t98 * t104, -t110, 0, 0, 0; 0, t104 * t117, t103 * t104 + t105 * t118, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (30->18), mult. (87->31), div. (0->0), fcn. (134->8), ass. (0->25)
	t94 = sin(pkin(6));
	t97 = sin(qJ(2));
	t111 = t94 * t97;
	t99 = cos(qJ(3));
	t110 = t94 * t99;
	t98 = sin(qJ(1));
	t109 = t98 * t97;
	t100 = cos(qJ(2));
	t108 = t100 * t94;
	t101 = cos(qJ(1));
	t107 = t101 * t94;
	t106 = t101 * t97;
	t105 = t98 * t100;
	t104 = t101 * t100;
	t95 = cos(pkin(6));
	t90 = t95 * t106 + t105;
	t96 = sin(qJ(3));
	t103 = t99 * t107 + t90 * t96;
	t102 = -t96 * t107 + t90 * t99;
	t92 = -t95 * t109 + t104;
	t91 = -t95 * t105 - t106;
	t89 = -t95 * t104 + t109;
	t88 = t98 * t94 * t96 + t92 * t99;
	t87 = -t98 * t110 + t92 * t96;
	t1 = [-t103, t91 * t96, t88, 0, 0, 0; t87, -t89 * t96, t102, 0, 0, 0; 0, t96 * t108, t97 * t110 + t95 * t96, 0, 0, 0; t102, -t91 * t99, t87, 0, 0, 0; -t88, t89 * t99, t103, 0, 0, 0; 0, -t99 * t108, t96 * t111 - t95 * t99, 0, 0, 0; t89, -t92, 0, 0, 0, 0; t91, -t90, 0, 0, 0, 0; 0, -t111, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (74->31), mult. (219->64), div. (0->0), fcn. (320->10), ass. (0->34)
	t134 = sin(qJ(2));
	t135 = sin(qJ(1));
	t138 = cos(qJ(2));
	t139 = cos(qJ(1));
	t152 = cos(pkin(6));
	t141 = t139 * t152;
	t126 = t134 * t141 + t135 * t138;
	t133 = sin(qJ(3));
	t137 = cos(qJ(3));
	t131 = sin(pkin(6));
	t147 = t131 * t139;
	t118 = t126 * t133 + t137 * t147;
	t125 = t135 * t134 - t138 * t141;
	t132 = sin(qJ(6));
	t136 = cos(qJ(6));
	t154 = t118 * t132 + t125 * t136;
	t153 = -t118 * t136 + t125 * t132;
	t149 = t131 * t133;
	t148 = t131 * t137;
	t146 = t132 * t133;
	t145 = t132 * t138;
	t144 = t133 * t136;
	t143 = t136 * t138;
	t142 = t135 * t152;
	t140 = -t126 * t137 + t133 * t147;
	t128 = -t134 * t142 + t139 * t138;
	t127 = -t139 * t134 - t138 * t142;
	t124 = t152 * t133 + t134 * t148;
	t123 = t134 * t149 - t152 * t137;
	t122 = t128 * t137 + t135 * t149;
	t121 = t128 * t133 - t135 * t148;
	t117 = t121 * t136 + t127 * t132;
	t116 = -t121 * t132 + t127 * t136;
	t1 = [t153, t127 * t144 - t128 * t132, t122 * t136, 0, 0, t116; t117, -t125 * t144 - t126 * t132, -t140 * t136, 0, 0, -t154; 0, (-t132 * t134 + t133 * t143) * t131, t124 * t136, 0, 0, -t123 * t132 + t131 * t143; t154, -t127 * t146 - t128 * t136, -t122 * t132, 0, 0, -t117; t116, t125 * t146 - t126 * t136, t140 * t132, 0, 0, t153; 0, (-t133 * t145 - t134 * t136) * t131, -t124 * t132, 0, 0, -t123 * t136 - t131 * t145; t140, t127 * t137, -t121, 0, 0, 0; t122, -t125 * t137, -t118, 0, 0, 0; 0, t138 * t148, -t123, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end