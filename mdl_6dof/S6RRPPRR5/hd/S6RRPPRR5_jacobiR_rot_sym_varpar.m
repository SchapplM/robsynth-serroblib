% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:43
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
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
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
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
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (8->6), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t61 = sin(qJ(2));
	t62 = sin(qJ(1));
	t68 = t62 * t61;
	t63 = cos(qJ(2));
	t67 = t62 * t63;
	t64 = cos(qJ(1));
	t66 = t64 * t61;
	t65 = t64 * t63;
	t60 = cos(pkin(6));
	t59 = sin(pkin(6));
	t58 = -t60 * t68 + t65;
	t57 = t60 * t67 + t66;
	t56 = t60 * t66 + t67;
	t55 = t60 * t65 - t68;
	t1 = [-t56, -t57, 0, 0, 0, 0; t58, t55, 0, 0, 0, 0; 0, t59 * t63, 0, 0, 0, 0; t64 * t59, 0, 0, 0, 0, 0; t62 * t59, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t55, t58, 0, 0, 0, 0; t57, t56, 0, 0, 0, 0; 0, t59 * t61, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->8), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t55 = sin(qJ(2));
	t56 = sin(qJ(1));
	t62 = t56 * t55;
	t57 = cos(qJ(2));
	t61 = t56 * t57;
	t58 = cos(qJ(1));
	t60 = t58 * t55;
	t59 = t58 * t57;
	t54 = cos(pkin(6));
	t53 = sin(pkin(6));
	t52 = -t54 * t62 + t59;
	t51 = t54 * t61 + t60;
	t50 = t54 * t60 + t61;
	t49 = t54 * t59 - t62;
	t1 = [-t50, -t51, 0, 0, 0, 0; t52, t49, 0, 0, 0, 0; 0, t53 * t57, 0, 0, 0, 0; t49, t52, 0, 0, 0, 0; t51, t50, 0, 0, 0, 0; 0, t53 * t55, 0, 0, 0, 0; -t58 * t53, 0, 0, 0, 0, 0; -t56 * t53, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (30->18), mult. (87->31), div. (0->0), fcn. (134->8), ass. (0->25)
	t83 = sin(pkin(6));
	t86 = sin(qJ(2));
	t100 = t83 * t86;
	t88 = cos(qJ(5));
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
	t85 = sin(qJ(5));
	t92 = -t79 * t85 + t88 * t97;
	t91 = -t79 * t88 - t85 * t97;
	t81 = -t84 * t96 + t93;
	t80 = -t84 * t95 - t94;
	t78 = -t84 * t93 + t96;
	t77 = -t87 * t83 * t85 + t81 * t88;
	t76 = -t81 * t85 - t87 * t99;
	t1 = [t91, t80 * t88, 0, 0, t76, 0; t77, -t78 * t88, 0, 0, t92, 0; 0, t88 * t98, 0, 0, -t85 * t100 - t84 * t88, 0; -t92, -t80 * t85, 0, 0, -t77, 0; t76, t78 * t85, 0, 0, t91, 0; 0, -t85 * t98, 0, 0, t84 * t85 - t86 * t99, 0; t78, -t81, 0, 0, 0, 0; t80, -t79, 0, 0, 0, 0; 0, -t100, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:07
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (74->28), mult. (219->64), div. (0->0), fcn. (320->10), ass. (0->34)
	t133 = sin(qJ(2));
	t134 = sin(qJ(1));
	t137 = cos(qJ(2));
	t138 = cos(qJ(1));
	t150 = cos(pkin(6));
	t139 = t138 * t150;
	t124 = t133 * t139 + t134 * t137;
	t132 = sin(qJ(5));
	t136 = cos(qJ(5));
	t130 = sin(pkin(6));
	t144 = t130 * t138;
	t117 = t124 * t136 + t132 * t144;
	t123 = t134 * t133 - t137 * t139;
	t131 = sin(qJ(6));
	t135 = cos(qJ(6));
	t152 = t117 * t131 + t123 * t135;
	t151 = -t117 * t135 + t123 * t131;
	t147 = t130 * t132;
	t146 = t130 * t136;
	t145 = t130 * t137;
	t143 = t131 * t136;
	t142 = t135 * t136;
	t141 = t136 * t137;
	t140 = t134 * t150;
	t116 = -t124 * t132 + t136 * t144;
	t126 = -t133 * t140 + t138 * t137;
	t125 = -t138 * t133 - t137 * t140;
	t122 = -t150 * t132 + t133 * t146;
	t121 = -t133 * t147 - t150 * t136;
	t120 = t126 * t136 - t134 * t147;
	t119 = t126 * t132 + t134 * t146;
	t115 = t120 * t135 + t125 * t131;
	t114 = -t120 * t131 + t125 * t135;
	t1 = [t151, t125 * t142 - t126 * t131, 0, 0, -t119 * t135, t114; t115, -t123 * t142 - t124 * t131, 0, 0, t116 * t135, -t152; 0, (-t131 * t133 + t135 * t141) * t130, 0, 0, t121 * t135, -t122 * t131 + t135 * t145; t152, -t125 * t143 - t126 * t135, 0, 0, t119 * t131, -t115; t114, t123 * t143 - t124 * t135, 0, 0, -t116 * t131, t151; 0, (-t131 * t141 - t133 * t135) * t130, 0, 0, -t121 * t131, -t122 * t135 - t131 * t145; t116, t125 * t132, 0, 0, t120, 0; t119, -t123 * t132, 0, 0, t117, 0; 0, t132 * t145, 0, 0, t122, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end