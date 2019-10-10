% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:16
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR14_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR14_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:16:27
	% EndTime: 2019-10-10 12:16:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:16:27
	% EndTime: 2019-10-10 12:16:27
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
	% StartTime: 2019-10-10 12:16:27
	% EndTime: 2019-10-10 12:16:27
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
	% StartTime: 2019-10-10 12:16:27
	% EndTime: 2019-10-10 12:16:27
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
	% StartTime: 2019-10-10 12:16:28
	% EndTime: 2019-10-10 12:16:28
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
	% StartTime: 2019-10-10 12:16:28
	% EndTime: 2019-10-10 12:16:28
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (71->31), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->36)
	t126 = cos(pkin(6));
	t129 = sin(qJ(2));
	t134 = cos(qJ(1));
	t137 = t134 * t129;
	t130 = sin(qJ(1));
	t133 = cos(qJ(2));
	t139 = t130 * t133;
	t121 = t126 * t137 + t139;
	t128 = sin(qJ(3));
	t132 = cos(qJ(3));
	t125 = sin(pkin(6));
	t144 = t125 * t134;
	t113 = t121 * t128 + t132 * t144;
	t136 = t134 * t133;
	t140 = t130 * t129;
	t120 = -t126 * t136 + t140;
	t127 = sin(qJ(5));
	t131 = cos(qJ(5));
	t150 = -t113 * t127 - t120 * t131;
	t149 = t113 * t131 - t120 * t127;
	t146 = t125 * t128;
	t145 = t125 * t132;
	t143 = t127 * t128;
	t142 = t127 * t133;
	t141 = t128 * t131;
	t138 = t131 * t133;
	t135 = -t121 * t132 + t128 * t144;
	t123 = -t126 * t140 + t136;
	t122 = t126 * t139 + t137;
	t119 = t126 * t128 + t129 * t145;
	t118 = -t126 * t132 + t129 * t146;
	t117 = t123 * t132 + t130 * t146;
	t116 = t123 * t128 - t130 * t145;
	t112 = t116 * t127 + t122 * t131;
	t111 = t116 * t131 - t122 * t127;
	t1 = [t150, -t122 * t143 + t123 * t131, t117 * t127, 0, t111, 0; t112, -t120 * t143 + t121 * t131, -t135 * t127, 0, t149, 0; 0, (t128 * t142 + t129 * t131) * t125, t119 * t127, 0, t118 * t131 + t125 * t142, 0; -t149, -t122 * t141 - t123 * t127, t117 * t131, 0, -t112, 0; t111, -t120 * t141 - t121 * t127, -t135 * t131, 0, t150, 0; 0, (-t127 * t129 + t128 * t138) * t125, t119 * t131, 0, -t118 * t127 + t125 * t138, 0; t135, -t122 * t132, -t116, 0, 0, 0; t117, -t120 * t132, -t113, 0, 0, 0; 0, t133 * t145, -t118, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:16:28
	% EndTime: 2019-10-10 12:16:28
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (141->33), mult. (275->62), div. (0->0), fcn. (402->10), ass. (0->39)
	t153 = cos(pkin(6));
	t155 = sin(qJ(2));
	t159 = cos(qJ(1));
	t162 = t159 * t155;
	t156 = sin(qJ(1));
	t158 = cos(qJ(2));
	t163 = t156 * t158;
	t145 = t153 * t162 + t163;
	t154 = sin(qJ(3));
	t157 = cos(qJ(3));
	t152 = sin(pkin(6));
	t166 = t152 * t159;
	t137 = t145 * t154 + t157 * t166;
	t161 = t159 * t158;
	t164 = t156 * t155;
	t144 = -t153 * t161 + t164;
	t151 = qJ(5) + qJ(6);
	t149 = sin(t151);
	t150 = cos(t151);
	t132 = -t137 * t149 - t144 * t150;
	t131 = t137 * t150 - t144 * t149;
	t171 = t149 * t154;
	t170 = t150 * t154;
	t169 = t152 * t154;
	t168 = t152 * t157;
	t167 = t152 * t158;
	t165 = t154 * t158;
	t160 = -t145 * t157 + t154 * t166;
	t147 = -t153 * t164 + t161;
	t146 = t153 * t163 + t162;
	t143 = t153 * t154 + t155 * t168;
	t142 = -t153 * t157 + t155 * t169;
	t141 = t147 * t157 + t156 * t169;
	t140 = t147 * t154 - t156 * t168;
	t136 = -t142 * t149 + t150 * t167;
	t135 = t142 * t150 + t149 * t167;
	t134 = t140 * t149 + t146 * t150;
	t133 = t140 * t150 - t146 * t149;
	t1 = [t132, -t146 * t171 + t147 * t150, t141 * t149, 0, t133, t133; t134, -t144 * t171 + t145 * t150, -t160 * t149, 0, t131, t131; 0, (t149 * t165 + t150 * t155) * t152, t143 * t149, 0, t135, t135; -t131, -t146 * t170 - t147 * t149, t141 * t150, 0, -t134, -t134; t133, -t144 * t170 - t145 * t149, -t160 * t150, 0, t132, t132; 0, (-t149 * t155 + t150 * t165) * t152, t143 * t150, 0, t136, t136; t160, -t146 * t157, -t140, 0, 0, 0; t141, -t144 * t157, -t137, 0, 0, 0; 0, t157 * t167, -t142, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end