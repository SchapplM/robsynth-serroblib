% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR5
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
% Datum: 2019-10-10 10:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
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
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
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
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
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
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
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
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (114->27), mult. (317->61), div. (0->0), fcn. (452->12), ass. (0->32)
	t146 = sin(pkin(12));
	t155 = cos(qJ(4));
	t165 = t146 * t155;
	t148 = sin(pkin(6));
	t154 = sin(qJ(1));
	t164 = t148 * t154;
	t157 = cos(qJ(1));
	t163 = t148 * t157;
	t149 = cos(pkin(12));
	t162 = t149 * t155;
	t151 = cos(pkin(6));
	t147 = sin(pkin(11));
	t150 = cos(pkin(11));
	t153 = sin(qJ(2));
	t156 = cos(qJ(2));
	t159 = t156 * t147 + t153 * t150;
	t141 = t159 * t151;
	t142 = t153 * t147 - t156 * t150;
	t161 = t157 * t141 - t154 * t142;
	t160 = -t154 * t141 - t157 * t142;
	t152 = sin(qJ(4));
	t127 = -t152 * t161 - t155 * t163;
	t128 = t152 * t163 - t155 * t161;
	t158 = t142 * t151;
	t140 = t159 * t148;
	t139 = t142 * t148;
	t137 = -t140 * t152 + t151 * t155;
	t135 = t154 * t158 - t157 * t159;
	t132 = -t154 * t159 - t157 * t158;
	t130 = t152 * t164 + t155 * t160;
	t129 = t152 * t160 - t155 * t164;
	t1 = [t128 * t149 + t132 * t146, t135 * t162 + t146 * t160, 0, -t129 * t149, 0, 0; t130 * t149 - t135 * t146, t132 * t162 + t146 * t161, 0, t127 * t149, 0, 0; 0, -t139 * t162 + t140 * t146, 0, t137 * t149, 0, 0; -t128 * t146 + t132 * t149, -t135 * t165 + t149 * t160, 0, t129 * t146, 0, 0; -t130 * t146 - t135 * t149, -t132 * t165 + t149 * t161, 0, -t127 * t146, 0, 0; 0, t139 * t165 + t140 * t149, 0, -t137 * t146, 0, 0; t127, t135 * t152, 0, t130, 0, 0; t129, t132 * t152, 0, -t128, 0, 0; 0, -t139 * t152, 0, t140 * t155 + t151 * t152, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (192->32), mult. (425->65), div. (0->0), fcn. (606->12), ass. (0->38)
	t176 = sin(qJ(4));
	t179 = cos(qJ(4));
	t175 = cos(pkin(6));
	t172 = sin(pkin(11));
	t174 = cos(pkin(11));
	t177 = sin(qJ(2));
	t180 = cos(qJ(2));
	t183 = t180 * t172 + t177 * t174;
	t163 = t183 * t175;
	t164 = t177 * t172 - t180 * t174;
	t178 = sin(qJ(1));
	t181 = cos(qJ(1));
	t185 = t181 * t163 - t178 * t164;
	t173 = sin(pkin(6));
	t188 = t173 * t181;
	t149 = t176 * t188 - t179 * t185;
	t182 = t164 * t175;
	t153 = -t178 * t183 - t181 * t182;
	t171 = pkin(12) + qJ(6);
	t169 = sin(t171);
	t170 = cos(t171);
	t195 = t149 * t169 - t153 * t170;
	t194 = t149 * t170 + t153 * t169;
	t191 = t169 * t179;
	t190 = t170 * t179;
	t189 = t173 * t178;
	t184 = -t178 * t163 - t181 * t164;
	t147 = -t176 * t185 - t179 * t188;
	t162 = t183 * t173;
	t161 = t164 * t173;
	t159 = t162 * t179 + t175 * t176;
	t158 = -t162 * t176 + t175 * t179;
	t156 = t178 * t182 - t181 * t183;
	t151 = t176 * t189 + t179 * t184;
	t150 = t176 * t184 - t179 * t189;
	t146 = t151 * t170 - t156 * t169;
	t145 = -t151 * t169 - t156 * t170;
	t1 = [t194, t156 * t190 + t169 * t184, 0, -t150 * t170, 0, t145; t146, t153 * t190 + t169 * t185, 0, t147 * t170, 0, t195; 0, -t161 * t190 + t162 * t169, 0, t158 * t170, 0, -t159 * t169 + t161 * t170; -t195, -t156 * t191 + t170 * t184, 0, t150 * t169, 0, -t146; t145, -t153 * t191 + t170 * t185, 0, -t147 * t169, 0, t194; 0, t161 * t191 + t162 * t170, 0, -t158 * t169, 0, -t159 * t170 - t161 * t169; t147, t156 * t176, 0, t151, 0, 0; t150, t153 * t176, 0, -t149, 0, 0; 0, -t161 * t176, 0, t159, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end