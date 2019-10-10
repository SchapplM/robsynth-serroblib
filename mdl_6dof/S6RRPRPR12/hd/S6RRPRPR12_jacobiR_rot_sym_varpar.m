% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR12_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR12_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
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
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
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
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
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
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (52->16), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->26)
	t95 = sin(pkin(6));
	t97 = sin(qJ(2));
	t110 = t95 * t97;
	t98 = sin(qJ(1));
	t109 = t95 * t98;
	t99 = cos(qJ(2));
	t108 = t95 * t99;
	t107 = t98 * t97;
	t106 = t98 * t99;
	t100 = cos(qJ(1));
	t105 = t100 * t95;
	t104 = t100 * t97;
	t103 = t100 * t99;
	t96 = cos(pkin(6));
	t86 = -t96 * t103 + t107;
	t94 = qJ(4) + pkin(11);
	t92 = sin(t94);
	t93 = cos(t94);
	t102 = t93 * t105 - t86 * t92;
	t101 = t92 * t105 + t86 * t93;
	t89 = -t96 * t107 + t103;
	t88 = t96 * t106 + t104;
	t87 = t96 * t104 + t106;
	t85 = t93 * t109 + t88 * t92;
	t84 = -t92 * t109 + t88 * t93;
	t1 = [t102, t89 * t92, 0, t84, 0, 0; t85, t87 * t92, 0, t101, 0, 0; 0, t92 * t110, 0, -t93 * t108 - t96 * t92, 0, 0; -t101, t89 * t93, 0, -t85, 0, 0; t84, t87 * t93, 0, t102, 0, 0; 0, t93 * t110, 0, t92 * t108 - t96 * t93, 0, 0; -t87, -t88, 0, 0, 0, 0; t89, -t86, 0, 0, 0, 0; 0, t108, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (128->32), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
	t142 = cos(pkin(6));
	t147 = cos(qJ(2));
	t148 = cos(qJ(1));
	t149 = t148 * t147;
	t144 = sin(qJ(2));
	t145 = sin(qJ(1));
	t152 = t145 * t144;
	t132 = -t142 * t149 + t152;
	t140 = qJ(4) + pkin(11);
	t138 = sin(t140);
	t139 = cos(t140);
	t141 = sin(pkin(6));
	t155 = t141 * t148;
	t127 = -t132 * t138 + t139 * t155;
	t150 = t148 * t144;
	t151 = t145 * t147;
	t133 = t142 * t150 + t151;
	t143 = sin(qJ(6));
	t146 = cos(qJ(6));
	t163 = t127 * t143 + t133 * t146;
	t162 = t127 * t146 - t133 * t143;
	t159 = t138 * t143;
	t158 = t138 * t146;
	t157 = t141 * t145;
	t156 = t141 * t147;
	t154 = t143 * t144;
	t153 = t144 * t146;
	t126 = t132 * t139 + t138 * t155;
	t135 = -t142 * t152 + t149;
	t134 = t142 * t151 + t150;
	t131 = -t138 * t156 + t142 * t139;
	t130 = -t142 * t138 - t139 * t156;
	t125 = t134 * t138 + t139 * t157;
	t124 = -t134 * t139 + t138 * t157;
	t123 = t125 * t146 + t135 * t143;
	t122 = -t125 * t143 + t135 * t146;
	t1 = [t162, -t134 * t143 + t135 * t158, 0, -t124 * t146, 0, t122; t123, -t132 * t143 + t133 * t158, 0, t126 * t146, 0, t163; 0, (t138 * t153 + t143 * t147) * t141, 0, t130 * t146, 0, -t131 * t143 + t141 * t153; -t163, -t134 * t146 - t135 * t159, 0, t124 * t143, 0, -t123; t122, -t132 * t146 - t133 * t159, 0, -t126 * t143, 0, t162; 0, (-t138 * t154 + t146 * t147) * t141, 0, -t130 * t143, 0, -t131 * t146 - t141 * t154; t126, -t135 * t139, 0, t125, 0, 0; t124, -t133 * t139, 0, -t127, 0, 0; 0, -t141 * t144 * t139, 0, t131, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end