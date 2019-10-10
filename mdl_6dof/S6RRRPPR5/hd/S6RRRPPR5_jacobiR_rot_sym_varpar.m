% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
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
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
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
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
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
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (55->16), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->26)
	t91 = sin(pkin(6));
	t93 = sin(qJ(2));
	t106 = t91 * t93;
	t94 = sin(qJ(1));
	t105 = t91 * t94;
	t95 = cos(qJ(2));
	t104 = t91 * t95;
	t96 = cos(qJ(1));
	t103 = t91 * t96;
	t102 = t94 * t93;
	t101 = t94 * t95;
	t100 = t96 * t93;
	t99 = t96 * t95;
	t92 = cos(pkin(6));
	t84 = t92 * t100 + t101;
	t90 = qJ(3) + pkin(11);
	t88 = sin(t90);
	t89 = cos(t90);
	t98 = t88 * t103 - t84 * t89;
	t97 = t89 * t103 + t84 * t88;
	t86 = -t92 * t102 + t99;
	t85 = t92 * t101 + t100;
	t83 = t92 * t99 - t102;
	t82 = t88 * t105 + t86 * t89;
	t81 = t89 * t105 - t86 * t88;
	t1 = [t98, -t85 * t89, t81, 0, 0, 0; t82, t83 * t89, -t97, 0, 0, 0; 0, t89 * t104, -t88 * t106 + t92 * t89, 0, 0, 0; t97, t85 * t88, -t82, 0, 0, 0; t81, -t83 * t88, t98, 0, 0, 0; 0, -t88 * t104, -t89 * t106 - t92 * t88, 0, 0, 0; t83, t86, 0, 0, 0, 0; t85, t84, 0, 0, 0, 0; 0, t106, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (93->26), mult. (163->58), div. (0->0), fcn. (238->10), ass. (0->31)
	t121 = qJ(3) + pkin(11);
	t120 = cos(t121);
	t122 = sin(pkin(12));
	t139 = t120 * t122;
	t124 = cos(pkin(12));
	t138 = t120 * t124;
	t128 = cos(qJ(2));
	t137 = t120 * t128;
	t123 = sin(pkin(6));
	t126 = sin(qJ(2));
	t136 = t123 * t126;
	t127 = sin(qJ(1));
	t135 = t123 * t127;
	t129 = cos(qJ(1));
	t134 = t123 * t129;
	t133 = t127 * t126;
	t132 = t127 * t128;
	t131 = t129 * t126;
	t130 = t129 * t128;
	t125 = cos(pkin(6));
	t115 = t125 * t131 + t132;
	t119 = sin(t121);
	t109 = -t115 * t119 - t120 * t134;
	t110 = -t115 * t120 + t119 * t134;
	t117 = -t125 * t133 + t130;
	t116 = t125 * t132 + t131;
	t114 = t125 * t130 - t133;
	t113 = -t119 * t136 + t125 * t120;
	t112 = t117 * t120 + t119 * t135;
	t111 = t117 * t119 - t120 * t135;
	t1 = [t110 * t124 + t114 * t122, -t116 * t138 + t117 * t122, -t111 * t124, 0, 0, 0; t112 * t124 + t116 * t122, t114 * t138 + t115 * t122, t109 * t124, 0, 0, 0; 0, (t122 * t126 + t124 * t137) * t123, t113 * t124, 0, 0, 0; -t110 * t122 + t114 * t124, t116 * t139 + t117 * t124, t111 * t122, 0, 0, 0; -t112 * t122 + t116 * t124, -t114 * t139 + t115 * t124, -t109 * t122, 0, 0, 0; 0, (-t122 * t137 + t124 * t126) * t123, -t113 * t122, 0, 0, 0; t109, -t116 * t119, t112, 0, 0, 0; t111, t114 * t119, -t110, 0, 0, 0; 0, t123 * t128 * t119, t125 * t119 + t120 * t136, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:24:03
	% EndTime: 2019-10-10 11:24:03
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (163->32), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->38)
	t147 = cos(pkin(6));
	t148 = sin(qJ(2));
	t151 = cos(qJ(1));
	t153 = t151 * t148;
	t149 = sin(qJ(1));
	t150 = cos(qJ(2));
	t154 = t149 * t150;
	t135 = t147 * t153 + t154;
	t145 = qJ(3) + pkin(11);
	t141 = sin(t145);
	t143 = cos(t145);
	t146 = sin(pkin(6));
	t156 = t146 * t151;
	t129 = -t135 * t143 + t141 * t156;
	t152 = t151 * t150;
	t155 = t149 * t148;
	t134 = -t147 * t152 + t155;
	t144 = pkin(12) + qJ(6);
	t140 = sin(t144);
	t142 = cos(t144);
	t166 = t129 * t140 + t134 * t142;
	t165 = t129 * t142 - t134 * t140;
	t162 = t140 * t143;
	t161 = t142 * t143;
	t160 = t143 * t150;
	t159 = t146 * t148;
	t158 = t146 * t149;
	t157 = t146 * t150;
	t127 = -t135 * t141 - t143 * t156;
	t137 = -t147 * t155 + t152;
	t136 = t147 * t154 + t153;
	t133 = t147 * t141 + t143 * t159;
	t132 = -t141 * t159 + t147 * t143;
	t131 = t137 * t143 + t141 * t158;
	t130 = t137 * t141 - t143 * t158;
	t126 = t131 * t142 + t136 * t140;
	t125 = -t131 * t140 + t136 * t142;
	t1 = [t165, -t136 * t161 + t137 * t140, -t130 * t142, 0, 0, t125; t126, -t134 * t161 + t135 * t140, t127 * t142, 0, 0, t166; 0, (t140 * t148 + t142 * t160) * t146, t132 * t142, 0, 0, -t133 * t140 - t142 * t157; -t166, t136 * t162 + t137 * t142, t130 * t140, 0, 0, -t126; t125, t134 * t162 + t135 * t142, -t127 * t140, 0, 0, t165; 0, (-t140 * t160 + t142 * t148) * t146, -t132 * t140, 0, 0, -t133 * t142 + t140 * t157; t127, -t136 * t141, t131, 0, 0, 0; t130, -t134 * t141, -t129, 0, 0, 0; 0, t141 * t157, t133, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end