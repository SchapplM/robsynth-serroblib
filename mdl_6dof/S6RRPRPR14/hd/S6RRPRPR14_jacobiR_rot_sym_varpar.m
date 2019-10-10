% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:28
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR14_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR14_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
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
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:15
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
	% StartTime: 2019-10-10 10:28:15
	% EndTime: 2019-10-10 10:28:16
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t62 = sin(qJ(2));
	t63 = sin(qJ(1));
	t69 = t63 * t62;
	t64 = cos(qJ(2));
	t68 = t63 * t64;
	t65 = cos(qJ(1));
	t67 = t65 * t62;
	t66 = t65 * t64;
	t61 = cos(pkin(6));
	t60 = sin(pkin(6));
	t59 = -t61 * t69 + t66;
	t58 = t61 * t68 + t67;
	t57 = t61 * t67 + t68;
	t56 = -t61 * t66 + t69;
	t1 = [t65 * t60, 0, 0, 0, 0, 0; t63 * t60, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t57, t58, 0, 0, 0, 0; -t59, t56, 0, 0, 0, 0; 0, -t60 * t64, 0, 0, 0, 0; -t56, t59, 0, 0, 0, 0; t58, t57, 0, 0, 0, 0; 0, t60 * t62, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:16
	% EndTime: 2019-10-10 10:28:16
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (26->15), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->25)
	t85 = sin(pkin(6));
	t87 = sin(qJ(4));
	t102 = t85 * t87;
	t90 = cos(qJ(4));
	t101 = t85 * t90;
	t91 = cos(qJ(2));
	t100 = t85 * t91;
	t92 = cos(qJ(1));
	t99 = t85 * t92;
	t88 = sin(qJ(2));
	t89 = sin(qJ(1));
	t98 = t89 * t88;
	t97 = t89 * t91;
	t96 = t92 * t88;
	t95 = t92 * t91;
	t86 = cos(pkin(6));
	t79 = -t86 * t95 + t98;
	t94 = -t79 * t87 + t90 * t99;
	t93 = t79 * t90 + t87 * t99;
	t82 = -t86 * t98 + t95;
	t81 = t86 * t97 + t96;
	t80 = t86 * t96 + t97;
	t78 = t101 * t89 + t81 * t87;
	t77 = -t102 * t89 + t81 * t90;
	t1 = [t94, t82 * t87, 0, t77, 0, 0; t78, t80 * t87, 0, t93, 0, 0; 0, t88 * t102, 0, -t100 * t90 - t86 * t87, 0, 0; -t93, t82 * t90, 0, -t78, 0, 0; t77, t80 * t90, 0, t94, 0, 0; 0, t88 * t101, 0, t100 * t87 - t86 * t90, 0, 0; -t80, -t81, 0, 0, 0, 0; t82, -t79, 0, 0, 0, 0; 0, t100, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:16
	% EndTime: 2019-10-10 10:28:16
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (32->21), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->25)
	t100 = sin(pkin(6));
	t102 = sin(qJ(4));
	t117 = t100 * t102;
	t105 = cos(qJ(4));
	t116 = t100 * t105;
	t106 = cos(qJ(2));
	t115 = t100 * t106;
	t107 = cos(qJ(1));
	t114 = t100 * t107;
	t103 = sin(qJ(2));
	t104 = sin(qJ(1));
	t113 = t104 * t103;
	t112 = t104 * t106;
	t111 = t107 * t103;
	t110 = t107 * t106;
	t101 = cos(pkin(6));
	t95 = -t101 * t110 + t113;
	t109 = t95 * t102 - t105 * t114;
	t108 = t102 * t114 + t95 * t105;
	t98 = -t101 * t113 + t110;
	t97 = t101 * t112 + t111;
	t96 = t101 * t111 + t112;
	t94 = t97 * t102 + t104 * t116;
	t93 = t104 * t117 - t97 * t105;
	t1 = [-t96, -t97, 0, 0, 0, 0; t98, -t95, 0, 0, 0, 0; 0, t115, 0, 0, 0, 0; t109, -t98 * t102, 0, t93, 0, 0; -t94, -t96 * t102, 0, -t108, 0, 0; 0, -t103 * t117, 0, t101 * t102 + t105 * t115, 0, 0; t108, -t98 * t105, 0, t94, 0, 0; t93, -t96 * t105, 0, t109, 0, 0; 0, -t103 * t116, 0, t101 * t105 - t102 * t115, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:28:16
	% EndTime: 2019-10-10 10:28:16
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (71->28), mult. (219->61), div. (0->0), fcn. (320->10), ass. (0->35)
	t130 = cos(pkin(6));
	t133 = sin(qJ(2));
	t138 = cos(qJ(1));
	t142 = t138 * t133;
	t134 = sin(qJ(1));
	t137 = cos(qJ(2));
	t144 = t134 * t137;
	t125 = t130 * t142 + t144;
	t131 = sin(qJ(6));
	t135 = cos(qJ(6));
	t141 = t138 * t137;
	t145 = t134 * t133;
	t124 = -t130 * t141 + t145;
	t132 = sin(qJ(4));
	t136 = cos(qJ(4));
	t129 = sin(pkin(6));
	t147 = t129 * t138;
	t139 = t124 * t136 + t132 * t147;
	t154 = -t125 * t135 + t131 * t139;
	t153 = t125 * t131 + t135 * t139;
	t150 = t129 * t133;
	t149 = t129 * t134;
	t148 = t129 * t137;
	t146 = t131 * t136;
	t143 = t135 * t136;
	t140 = -t124 * t132 + t136 * t147;
	t127 = -t130 * t145 + t141;
	t126 = t130 * t144 + t142;
	t123 = t130 * t136 - t132 * t148;
	t122 = t130 * t132 + t136 * t148;
	t118 = t126 * t132 + t136 * t149;
	t117 = -t126 * t136 + t132 * t149;
	t116 = t117 * t131 + t127 * t135;
	t115 = t117 * t135 - t127 * t131;
	t1 = [t154, -t126 * t135 - t127 * t146, 0, t118 * t131, 0, t115; t116, -t124 * t135 - t125 * t146, 0, -t140 * t131, 0, -t153; 0, (-t133 * t146 + t135 * t137) * t129, 0, t123 * t131, 0, t122 * t135 - t131 * t150; t153, t126 * t131 - t127 * t143, 0, t118 * t135, 0, -t116; t115, t124 * t131 - t125 * t143, 0, -t140 * t135, 0, t154; 0, (-t131 * t137 - t133 * t143) * t129, 0, t123 * t135, 0, -t122 * t131 - t135 * t150; t140, t127 * t132, 0, -t117, 0, 0; t118, t125 * t132, 0, t139, 0, 0; 0, t132 * t150, 0, -t122, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end