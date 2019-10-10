% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
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
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
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
	% StartTime: 2019-10-10 10:13:32
	% EndTime: 2019-10-10 10:13:32
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (26->10), mult. (74->18), div. (0->0), fcn. (112->8), ass. (0->17)
	t70 = cos(pkin(6));
	t67 = sin(pkin(11));
	t69 = cos(pkin(11));
	t71 = sin(qJ(2));
	t73 = cos(qJ(2));
	t75 = t73 * t67 + t71 * t69;
	t63 = t75 * t70;
	t64 = t71 * t67 - t73 * t69;
	t72 = sin(qJ(1));
	t74 = cos(qJ(1));
	t77 = -t74 * t63 + t72 * t64;
	t76 = t72 * t63 + t74 * t64;
	t68 = sin(pkin(6));
	t62 = t64 * t70;
	t61 = t72 * t62 - t74 * t75;
	t60 = -t74 * t62 - t72 * t75;
	t1 = [t77, t61, 0, 0, 0, 0; -t76, t60, 0, 0, 0, 0; 0, -t64 * t68, 0, 0, 0, 0; -t60, t76, 0, 0, 0, 0; t61, t77, 0, 0, 0, 0; 0, -t75 * t68, 0, 0, 0, 0; t74 * t68, 0, 0, 0, 0, 0; t72 * t68, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:33
	% EndTime: 2019-10-10 10:13:33
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (66->18), mult. (182->36), div. (0->0), fcn. (266->10), ass. (0->27)
	t115 = sin(pkin(6));
	t120 = sin(qJ(1));
	t128 = t115 * t120;
	t123 = cos(qJ(1));
	t127 = t115 * t123;
	t117 = cos(pkin(6));
	t114 = sin(pkin(11));
	t116 = cos(pkin(11));
	t119 = sin(qJ(2));
	t122 = cos(qJ(2));
	t125 = t122 * t114 + t119 * t116;
	t109 = t125 * t117;
	t110 = t119 * t114 - t122 * t116;
	t103 = t123 * t109 - t120 * t110;
	t118 = sin(qJ(4));
	t121 = cos(qJ(4));
	t126 = -t103 * t121 + t118 * t127;
	t105 = -t120 * t109 - t123 * t110;
	t124 = t103 * t118 + t121 * t127;
	t108 = t110 * t117;
	t107 = t125 * t115;
	t106 = t110 * t115;
	t104 = t120 * t108 - t123 * t125;
	t102 = -t123 * t108 - t120 * t125;
	t101 = t105 * t121 + t118 * t128;
	t100 = -t105 * t118 + t121 * t128;
	t1 = [t126, t104 * t121, 0, t100, 0, 0; t101, t102 * t121, 0, -t124, 0, 0; 0, -t106 * t121, 0, -t107 * t118 + t117 * t121, 0, 0; t124, -t104 * t118, 0, -t101, 0, 0; t100, -t102 * t118, 0, t126, 0, 0; 0, t106 * t118, 0, -t107 * t121 - t117 * t118, 0, 0; t102, t105, 0, 0, 0, 0; -t104, t103, 0, 0, 0, 0; 0, t107, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:33
	% EndTime: 2019-10-10 10:13:33
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (66->18), mult. (182->36), div. (0->0), fcn. (266->10), ass. (0->27)
	t127 = sin(pkin(6));
	t132 = sin(qJ(1));
	t140 = t127 * t132;
	t135 = cos(qJ(1));
	t139 = t127 * t135;
	t129 = cos(pkin(6));
	t126 = sin(pkin(11));
	t128 = cos(pkin(11));
	t131 = sin(qJ(2));
	t134 = cos(qJ(2));
	t138 = t134 * t126 + t131 * t128;
	t122 = t138 * t129;
	t123 = t131 * t126 - t134 * t128;
	t116 = t135 * t122 - t132 * t123;
	t118 = -t132 * t122 - t135 * t123;
	t130 = sin(qJ(4));
	t133 = cos(qJ(4));
	t137 = t116 * t130 + t133 * t139;
	t136 = t116 * t133 - t130 * t139;
	t121 = t123 * t129;
	t120 = t138 * t127;
	t119 = t123 * t127;
	t117 = t132 * t121 - t135 * t138;
	t115 = -t135 * t121 - t132 * t138;
	t114 = t118 * t133 + t130 * t140;
	t113 = t118 * t130 - t133 * t140;
	t1 = [t115, t118, 0, 0, 0, 0; -t117, t116, 0, 0, 0, 0; 0, t120, 0, 0, 0, 0; t136, -t117 * t133, 0, t113, 0, 0; -t114, -t115 * t133, 0, t137, 0, 0; 0, t119 * t133, 0, t120 * t130 - t129 * t133, 0, 0; -t137, t117 * t130, 0, t114, 0, 0; t113, t115 * t130, 0, t136, 0, 0; 0, -t119 * t130, 0, t120 * t133 + t129 * t130, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:13:33
	% EndTime: 2019-10-10 10:13:33
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (151->32), mult. (425->65), div. (0->0), fcn. (606->12), ass. (0->37)
	t167 = sin(qJ(4));
	t171 = cos(qJ(4));
	t165 = cos(pkin(6));
	t162 = sin(pkin(11));
	t164 = cos(pkin(11));
	t168 = sin(qJ(2));
	t172 = cos(qJ(2));
	t176 = t172 * t162 + t168 * t164;
	t157 = t176 * t165;
	t158 = t168 * t162 - t172 * t164;
	t169 = sin(qJ(1));
	t173 = cos(qJ(1));
	t178 = t173 * t157 - t169 * t158;
	t163 = sin(pkin(6));
	t183 = t163 * t173;
	t141 = t167 * t178 + t171 * t183;
	t174 = t158 * t165;
	t147 = -t169 * t176 - t173 * t174;
	t166 = sin(qJ(6));
	t170 = cos(qJ(6));
	t188 = -t141 * t166 + t147 * t170;
	t187 = t141 * t170 + t147 * t166;
	t184 = t163 * t169;
	t182 = t166 * t167;
	t181 = t167 * t170;
	t177 = -t169 * t157 - t173 * t158;
	t175 = t167 * t183 - t171 * t178;
	t156 = t176 * t163;
	t155 = t158 * t163;
	t153 = t156 * t171 + t165 * t167;
	t152 = t156 * t167 - t165 * t171;
	t150 = t169 * t174 - t173 * t176;
	t145 = t167 * t184 + t171 * t177;
	t144 = t167 * t177 - t171 * t184;
	t140 = t144 * t166 - t150 * t170;
	t139 = t144 * t170 + t150 * t166;
	t1 = [t188, t150 * t182 + t170 * t177, 0, t145 * t166, 0, t139; t140, t147 * t182 + t170 * t178, 0, -t175 * t166, 0, t187; 0, -t155 * t182 + t156 * t170, 0, t153 * t166, 0, t152 * t170 - t155 * t166; -t187, t150 * t181 - t166 * t177, 0, t145 * t170, 0, -t140; t139, t147 * t181 - t166 * t178, 0, -t175 * t170, 0, t188; 0, -t155 * t181 - t156 * t166, 0, t153 * t170, 0, -t152 * t166 - t155 * t170; t175, t150 * t171, 0, -t144, 0, 0; t145, t147 * t171, 0, -t141, 0, 0; 0, -t155 * t171, 0, -t152, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end