% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:33
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR10_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR10_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
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
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
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
	% StartTime: 2019-10-10 11:33:20
	% EndTime: 2019-10-10 11:33:20
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
	% StartTime: 2019-10-10 11:33:20
	% EndTime: 2019-10-10 11:33:20
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (29->15), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->25)
	t101 = sin(qJ(1));
	t97 = sin(pkin(6));
	t114 = t101 * t97;
	t103 = cos(qJ(2));
	t113 = t103 * t97;
	t104 = cos(qJ(1));
	t112 = t104 * t97;
	t100 = sin(qJ(2));
	t111 = t97 * t100;
	t110 = t101 * t100;
	t109 = t101 * t103;
	t108 = t104 * t100;
	t107 = t104 * t103;
	t102 = cos(qJ(3));
	t98 = cos(pkin(6));
	t94 = t98 * t108 + t109;
	t99 = sin(qJ(3));
	t106 = t94 * t102 - t99 * t112;
	t105 = t102 * t112 + t94 * t99;
	t96 = -t98 * t110 + t107;
	t95 = t98 * t109 + t108;
	t93 = t98 * t107 - t110;
	t92 = t96 * t102 + t99 * t114;
	t91 = -t102 * t114 + t96 * t99;
	t1 = [t93, t96, 0, 0, 0, 0; t95, t94, 0, 0, 0, 0; 0, t111, 0, 0, 0, 0; t106, t95 * t102, t91, 0, 0, 0; -t92, -t93 * t102, t105, 0, 0, 0; 0, -t102 * t113, -t98 * t102 + t99 * t111, 0, 0, 0; -t105, -t95 * t99, t92, 0, 0, 0; t91, t93 * t99, t106, 0, 0, 0; 0, t99 * t113, t102 * t111 + t98 * t99, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:20
	% EndTime: 2019-10-10 11:33:20
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (51->24), mult. (163->57), div. (0->0), fcn. (238->10), ass. (0->30)
	t107 = sin(pkin(11));
	t111 = sin(qJ(3));
	t127 = t107 * t111;
	t108 = sin(pkin(6));
	t126 = t108 * t111;
	t114 = cos(qJ(3));
	t125 = t108 * t114;
	t116 = cos(qJ(1));
	t124 = t108 * t116;
	t109 = cos(pkin(11));
	t123 = t109 * t111;
	t115 = cos(qJ(2));
	t122 = t111 * t115;
	t112 = sin(qJ(2));
	t113 = sin(qJ(1));
	t121 = t113 * t112;
	t120 = t113 * t115;
	t119 = t116 * t112;
	t118 = t116 * t115;
	t110 = cos(pkin(6));
	t104 = t110 * t119 + t120;
	t99 = -t104 * t111 - t114 * t124;
	t117 = -t104 * t114 + t111 * t124;
	t106 = -t110 * t121 + t118;
	t105 = t110 * t120 + t119;
	t103 = t110 * t118 - t121;
	t102 = t110 * t111 + t112 * t125;
	t101 = t106 * t114 + t113 * t126;
	t100 = t106 * t111 - t113 * t125;
	t1 = [t103 * t109 + t99 * t107, -t105 * t127 + t106 * t109, t101 * t107, 0, 0, 0; t100 * t107 + t105 * t109, t103 * t127 + t104 * t109, -t117 * t107, 0, 0, 0; 0, (t107 * t122 + t109 * t112) * t108, t102 * t107, 0, 0, 0; -t103 * t107 + t99 * t109, -t105 * t123 - t106 * t107, t101 * t109, 0, 0, 0; t100 * t109 - t105 * t107, t103 * t123 - t104 * t107, -t117 * t109, 0, 0, 0; 0, (-t107 * t112 + t109 * t122) * t108, t102 * t109, 0, 0, 0; t117, -t105 * t114, -t100, 0, 0, 0; t101, t103 * t114, t99, 0, 0, 0; 0, t115 * t125, t110 * t114 - t112 * t126, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:20
	% EndTime: 2019-10-10 11:33:20
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (109->32), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->37)
	t134 = cos(pkin(6));
	t136 = sin(qJ(2));
	t140 = cos(qJ(1));
	t143 = t140 * t136;
	t137 = sin(qJ(1));
	t139 = cos(qJ(2));
	t144 = t137 * t139;
	t126 = t134 * t143 + t144;
	t135 = sin(qJ(3));
	t138 = cos(qJ(3));
	t133 = sin(pkin(6));
	t147 = t133 * t140;
	t118 = t126 * t135 + t138 * t147;
	t142 = t140 * t139;
	t145 = t137 * t136;
	t125 = -t134 * t142 + t145;
	t132 = pkin(11) + qJ(6);
	t130 = sin(t132);
	t131 = cos(t132);
	t156 = -t118 * t130 - t125 * t131;
	t155 = t118 * t131 - t125 * t130;
	t152 = t130 * t135;
	t151 = t131 * t135;
	t150 = t133 * t135;
	t149 = t133 * t138;
	t148 = t133 * t139;
	t146 = t135 * t139;
	t141 = -t126 * t138 + t135 * t147;
	t128 = -t134 * t145 + t142;
	t127 = t134 * t144 + t143;
	t124 = t134 * t135 + t136 * t149;
	t123 = -t134 * t138 + t136 * t150;
	t122 = t128 * t138 + t137 * t150;
	t121 = t128 * t135 - t137 * t149;
	t117 = t121 * t130 + t127 * t131;
	t116 = t121 * t131 - t127 * t130;
	t1 = [t156, -t127 * t152 + t128 * t131, t122 * t130, 0, 0, t116; t117, -t125 * t152 + t126 * t131, -t141 * t130, 0, 0, t155; 0, (t130 * t146 + t131 * t136) * t133, t124 * t130, 0, 0, t123 * t131 + t130 * t148; -t155, -t127 * t151 - t128 * t130, t122 * t131, 0, 0, -t117; t116, -t125 * t151 - t126 * t130, -t141 * t131, 0, 0, t156; 0, (-t130 * t136 + t131 * t146) * t133, t124 * t131, 0, 0, -t123 * t130 + t131 * t148; t141, -t127 * t138, -t121, 0, 0, 0; t122, -t125 * t138, -t118, 0, 0, 0; 0, t138 * t148, -t123, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end