% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR9
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
% Datum: 2019-10-10 09:50
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR9_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR9_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
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
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
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
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
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
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (8->6), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t60 = sin(qJ(2));
	t61 = sin(qJ(1));
	t67 = t61 * t60;
	t62 = cos(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(1));
	t65 = t63 * t60;
	t64 = t63 * t62;
	t59 = cos(pkin(6));
	t58 = sin(pkin(6));
	t57 = -t59 * t67 + t64;
	t56 = t59 * t66 + t65;
	t55 = t59 * t65 + t66;
	t54 = t59 * t64 - t67;
	t1 = [t63 * t58, 0, 0, 0, 0, 0; t61 * t58, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t54, t57, 0, 0, 0, 0; t56, t55, 0, 0, 0, 0; 0, t58 * t60, 0, 0, 0, 0; -t55, -t56, 0, 0, 0, 0; t57, t54, 0, 0, 0, 0; 0, t58 * t62, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (27->17), mult. (87->31), div. (0->0), fcn. (134->8), ass. (0->25)
	t86 = sin(pkin(6));
	t89 = sin(qJ(2));
	t103 = t86 * t89;
	t91 = cos(qJ(5));
	t102 = t86 * t91;
	t92 = cos(qJ(2));
	t101 = t86 * t92;
	t93 = cos(qJ(1));
	t100 = t86 * t93;
	t90 = sin(qJ(1));
	t99 = t90 * t89;
	t98 = t90 * t92;
	t97 = t93 * t89;
	t96 = t93 * t92;
	t87 = cos(pkin(6));
	t81 = t87 * t97 + t98;
	t88 = sin(qJ(5));
	t95 = t91 * t100 - t81 * t88;
	t94 = t88 * t100 + t81 * t91;
	t83 = -t87 * t99 + t96;
	t82 = -t87 * t98 - t97;
	t80 = -t87 * t96 + t99;
	t79 = t90 * t102 + t83 * t88;
	t78 = -t90 * t86 * t88 + t83 * t91;
	t1 = [t95, t82 * t88, 0, 0, t78, 0; t79, -t80 * t88, 0, 0, t94, 0; 0, t88 * t101, 0, 0, t89 * t102 - t87 * t88, 0; -t94, t82 * t91, 0, 0, -t79, 0; t78, -t80 * t91, 0, 0, t95, 0; 0, t91 * t101, 0, 0, -t88 * t103 - t87 * t91, 0; t80, -t83, 0, 0, 0, 0; t82, -t81, 0, 0, 0, 0; 0, -t103, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (77->30), mult. (219->64), div. (0->0), fcn. (320->10), ass. (0->34)
	t137 = sin(qJ(2));
	t138 = sin(qJ(1));
	t141 = cos(qJ(2));
	t142 = cos(qJ(1));
	t154 = cos(pkin(6));
	t143 = t142 * t154;
	t128 = t137 * t143 + t138 * t141;
	t136 = sin(qJ(5));
	t140 = cos(qJ(5));
	t134 = sin(pkin(6));
	t149 = t134 * t142;
	t122 = -t128 * t136 + t140 * t149;
	t127 = t138 * t137 - t141 * t143;
	t135 = sin(qJ(6));
	t139 = cos(qJ(6));
	t156 = t122 * t135 - t127 * t139;
	t155 = t122 * t139 + t127 * t135;
	t151 = t134 * t136;
	t150 = t134 * t140;
	t148 = t135 * t136;
	t147 = t135 * t141;
	t146 = t136 * t139;
	t145 = t139 * t141;
	t144 = t138 * t154;
	t121 = t128 * t140 + t136 * t149;
	t130 = -t137 * t144 + t142 * t141;
	t129 = -t142 * t137 - t141 * t144;
	t126 = t137 * t151 + t154 * t140;
	t125 = -t154 * t136 + t137 * t150;
	t120 = t130 * t136 + t138 * t150;
	t119 = -t130 * t140 + t138 * t151;
	t118 = t120 * t139 + t129 * t135;
	t117 = -t120 * t135 + t129 * t139;
	t1 = [t155, t129 * t146 - t130 * t135, 0, 0, -t119 * t139, t117; t118, -t127 * t146 - t128 * t135, 0, 0, t121 * t139, t156; 0, (-t135 * t137 + t136 * t145) * t134, 0, 0, t125 * t139, -t126 * t135 + t134 * t145; -t156, -t129 * t148 - t130 * t139, 0, 0, t119 * t135, -t118; t117, t127 * t148 - t128 * t139, 0, 0, -t121 * t135, t155; 0, (-t136 * t147 - t137 * t139) * t134, 0, 0, -t125 * t135, -t126 * t139 - t134 * t147; t121, -t129 * t140, 0, 0, t120, 0; t119, t127 * t140, 0, 0, -t122, 0; 0, -t141 * t150, 0, 0, t126, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end