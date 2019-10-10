% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:33
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPP9_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP9_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:10
	% EndTime: 2019-10-10 12:33:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:10
	% EndTime: 2019-10-10 12:33:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:10
	% EndTime: 2019-10-10 12:33:10
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
	% StartTime: 2019-10-10 12:33:11
	% EndTime: 2019-10-10 12:33:11
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
	% StartTime: 2019-10-10 12:33:11
	% EndTime: 2019-10-10 12:33:11
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (74->30), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->36)
	t131 = cos(pkin(6));
	t134 = sin(qJ(2));
	t139 = cos(qJ(1));
	t141 = t139 * t134;
	t135 = sin(qJ(1));
	t138 = cos(qJ(2));
	t144 = t135 * t138;
	t125 = t131 * t141 + t144;
	t133 = sin(qJ(3));
	t137 = cos(qJ(3));
	t130 = sin(pkin(6));
	t147 = t130 * t139;
	t119 = -t125 * t137 + t133 * t147;
	t140 = t139 * t138;
	t145 = t135 * t134;
	t124 = -t131 * t140 + t145;
	t132 = sin(qJ(4));
	t136 = cos(qJ(4));
	t154 = t119 * t132 + t124 * t136;
	t153 = t119 * t136 - t124 * t132;
	t150 = t130 * t133;
	t149 = t130 * t137;
	t148 = t130 * t138;
	t146 = t132 * t137;
	t143 = t136 * t137;
	t142 = t137 * t138;
	t117 = -t125 * t133 - t137 * t147;
	t127 = -t131 * t145 + t140;
	t126 = t131 * t144 + t141;
	t123 = t131 * t133 + t134 * t149;
	t122 = t131 * t137 - t134 * t150;
	t121 = t127 * t137 + t135 * t150;
	t120 = t127 * t133 - t135 * t149;
	t116 = t121 * t136 + t126 * t132;
	t115 = -t121 * t132 + t126 * t136;
	t1 = [t153, -t126 * t143 + t127 * t132, -t120 * t136, t115, 0, 0; t116, -t124 * t143 + t125 * t132, t117 * t136, t154, 0, 0; 0, (t132 * t134 + t136 * t142) * t130, t122 * t136, -t123 * t132 - t136 * t148, 0, 0; -t154, t126 * t146 + t127 * t136, t120 * t132, -t116, 0, 0; t115, t124 * t146 + t125 * t136, -t117 * t132, t153, 0, 0; 0, (-t132 * t142 + t134 * t136) * t130, -t122 * t132, -t123 * t136 + t132 * t148, 0, 0; t117, -t126 * t133, t121, 0, 0, 0; t120, -t124 * t133, -t119, 0, 0, 0; 0, t133 * t148, t123, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:11
	% EndTime: 2019-10-10 12:33:11
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (74->30), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->36)
	t161 = cos(pkin(6));
	t164 = sin(qJ(2));
	t169 = cos(qJ(1));
	t171 = t169 * t164;
	t165 = sin(qJ(1));
	t168 = cos(qJ(2));
	t174 = t165 * t168;
	t155 = t161 * t171 + t174;
	t163 = sin(qJ(3));
	t167 = cos(qJ(3));
	t160 = sin(pkin(6));
	t177 = t160 * t169;
	t149 = -t155 * t167 + t163 * t177;
	t170 = t169 * t168;
	t175 = t165 * t164;
	t154 = -t161 * t170 + t175;
	t162 = sin(qJ(4));
	t166 = cos(qJ(4));
	t184 = t149 * t162 + t154 * t166;
	t183 = -t149 * t166 + t154 * t162;
	t180 = t160 * t163;
	t179 = t160 * t167;
	t178 = t160 * t168;
	t176 = t162 * t167;
	t173 = t166 * t167;
	t172 = t167 * t168;
	t147 = -t155 * t163 - t167 * t177;
	t157 = -t161 * t175 + t170;
	t156 = t161 * t174 + t171;
	t153 = t161 * t163 + t164 * t179;
	t152 = t161 * t167 - t164 * t180;
	t151 = t157 * t167 + t165 * t180;
	t150 = t157 * t163 - t165 * t179;
	t146 = t151 * t166 + t156 * t162;
	t145 = t151 * t162 - t156 * t166;
	t1 = [t147, -t156 * t163, t151, 0, 0, 0; t150, -t154 * t163, -t149, 0, 0, 0; 0, t163 * t178, t153, 0, 0, 0; t183, t156 * t173 - t157 * t162, t150 * t166, t145, 0, 0; -t146, t154 * t173 - t155 * t162, -t147 * t166, -t184, 0, 0; 0, (-t162 * t164 - t166 * t172) * t160, -t152 * t166, t153 * t162 + t166 * t178, 0, 0; t184, -t156 * t176 - t157 * t166, -t150 * t162, t146, 0, 0; t145, -t154 * t176 - t155 * t166, t147 * t162, t183, 0, 0; 0, (t162 * t172 - t164 * t166) * t160, t152 * t162, t153 * t166 - t162 * t178, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:11
	% EndTime: 2019-10-10 12:33:11
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (71->29), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->36)
	t158 = cos(pkin(6));
	t161 = sin(qJ(2));
	t166 = cos(qJ(1));
	t168 = t166 * t161;
	t162 = sin(qJ(1));
	t165 = cos(qJ(2));
	t171 = t162 * t165;
	t152 = t158 * t168 + t171;
	t160 = sin(qJ(3));
	t164 = cos(qJ(3));
	t157 = sin(pkin(6));
	t174 = t157 * t166;
	t146 = -t152 * t164 + t160 * t174;
	t167 = t166 * t165;
	t172 = t162 * t161;
	t151 = -t158 * t167 + t172;
	t159 = sin(qJ(4));
	t163 = cos(qJ(4));
	t181 = t146 * t159 + t151 * t163;
	t180 = t146 * t163 - t151 * t159;
	t177 = t157 * t160;
	t176 = t157 * t164;
	t175 = t157 * t165;
	t173 = t159 * t164;
	t170 = t163 * t164;
	t169 = t164 * t165;
	t144 = -t152 * t160 - t164 * t174;
	t154 = -t158 * t172 + t167;
	t153 = t158 * t171 + t168;
	t150 = t158 * t160 + t161 * t176;
	t149 = t158 * t164 - t161 * t177;
	t148 = t154 * t164 + t162 * t177;
	t147 = t154 * t160 - t162 * t176;
	t143 = t148 * t163 + t153 * t159;
	t142 = t148 * t159 - t153 * t163;
	t1 = [t144, -t153 * t160, t148, 0, 0, 0; t147, -t151 * t160, -t146, 0, 0, 0; 0, t160 * t175, t150, 0, 0, 0; t181, -t153 * t173 - t154 * t163, -t147 * t159, t143, 0, 0; t142, -t151 * t173 - t152 * t163, t144 * t159, -t180, 0, 0; 0, (t159 * t169 - t161 * t163) * t157, t149 * t159, t150 * t163 - t159 * t175, 0, 0; t180, -t153 * t170 + t154 * t159, -t147 * t163, -t142, 0, 0; t143, -t151 * t170 + t152 * t159, t144 * t163, t181, 0, 0; 0, (t159 * t161 + t163 * t169) * t157, t149 * t163, -t150 * t159 - t163 * t175, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end