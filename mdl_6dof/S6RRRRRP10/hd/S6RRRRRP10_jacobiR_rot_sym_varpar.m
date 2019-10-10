% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP10_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP10_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:11:34
	% EndTime: 2019-10-10 13:11:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:11:34
	% EndTime: 2019-10-10 13:11:34
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
	% StartTime: 2019-10-10 13:11:34
	% EndTime: 2019-10-10 13:11:34
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
	% StartTime: 2019-10-10 13:11:34
	% EndTime: 2019-10-10 13:11:34
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-10 13:11:34
	% EndTime: 2019-10-10 13:11:34
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
	% StartTime: 2019-10-10 13:11:35
	% EndTime: 2019-10-10 13:11:35
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (144->32), mult. (275->62), div. (0->0), fcn. (402->10), ass. (0->39)
	t154 = cos(pkin(6));
	t156 = sin(qJ(2));
	t160 = cos(qJ(1));
	t162 = t160 * t156;
	t157 = sin(qJ(1));
	t159 = cos(qJ(2));
	t164 = t157 * t159;
	t145 = t154 * t162 + t164;
	t155 = sin(qJ(3));
	t158 = cos(qJ(3));
	t153 = sin(pkin(6));
	t166 = t153 * t160;
	t139 = -t145 * t158 + t155 * t166;
	t161 = t160 * t159;
	t165 = t157 * t156;
	t144 = -t154 * t161 + t165;
	t152 = qJ(4) + qJ(5);
	t150 = sin(t152);
	t151 = cos(t152);
	t131 = t139 * t150 + t144 * t151;
	t132 = t139 * t151 - t144 * t150;
	t171 = t150 * t158;
	t170 = t151 * t158;
	t169 = t153 * t155;
	t168 = t153 * t158;
	t167 = t153 * t159;
	t163 = t158 * t159;
	t137 = -t145 * t155 - t158 * t166;
	t147 = -t154 * t165 + t161;
	t146 = t154 * t164 + t162;
	t143 = t154 * t155 + t156 * t168;
	t142 = t154 * t158 - t156 * t169;
	t141 = t147 * t158 + t157 * t169;
	t140 = t147 * t155 - t157 * t168;
	t136 = -t143 * t151 + t150 * t167;
	t135 = -t143 * t150 - t151 * t167;
	t134 = t141 * t151 + t146 * t150;
	t133 = -t141 * t150 + t146 * t151;
	t1 = [t132, -t146 * t170 + t147 * t150, -t140 * t151, t133, t133, 0; t134, -t144 * t170 + t145 * t150, t137 * t151, t131, t131, 0; 0, (t150 * t156 + t151 * t163) * t153, t142 * t151, t135, t135, 0; -t131, t146 * t171 + t147 * t151, t140 * t150, -t134, -t134, 0; t133, t144 * t171 + t145 * t151, -t137 * t150, t132, t132, 0; 0, (-t150 * t163 + t151 * t156) * t153, -t142 * t150, t136, t136, 0; t137, -t146 * t155, t141, 0, 0, 0; t140, -t144 * t155, -t139, 0, 0, 0; 0, t155 * t167, t143, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:11:35
	% EndTime: 2019-10-10 13:11:35
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (141->32), mult. (275->62), div. (0->0), fcn. (402->10), ass. (0->39)
	t196 = cos(pkin(6));
	t198 = sin(qJ(2));
	t202 = cos(qJ(1));
	t204 = t202 * t198;
	t199 = sin(qJ(1));
	t201 = cos(qJ(2));
	t206 = t199 * t201;
	t187 = t196 * t204 + t206;
	t197 = sin(qJ(3));
	t200 = cos(qJ(3));
	t195 = sin(pkin(6));
	t208 = t195 * t202;
	t180 = -t187 * t200 + t197 * t208;
	t203 = t202 * t201;
	t207 = t199 * t198;
	t186 = -t196 * t203 + t207;
	t194 = qJ(4) + qJ(5);
	t192 = sin(t194);
	t193 = cos(t194);
	t172 = t180 * t192 + t186 * t193;
	t216 = t180 * t193 - t186 * t192;
	t213 = t192 * t200;
	t212 = t193 * t200;
	t211 = t195 * t197;
	t210 = t195 * t200;
	t209 = t195 * t201;
	t205 = t200 * t201;
	t178 = -t187 * t197 - t200 * t208;
	t189 = -t196 * t207 + t203;
	t188 = t196 * t206 + t204;
	t185 = t196 * t197 + t198 * t210;
	t184 = t196 * t200 - t198 * t211;
	t182 = t189 * t200 + t199 * t211;
	t181 = t189 * t197 - t199 * t210;
	t177 = t185 * t193 - t192 * t209;
	t176 = -t185 * t192 - t193 * t209;
	t175 = t182 * t193 + t188 * t192;
	t174 = t182 * t192 - t188 * t193;
	t1 = [t216, -t188 * t212 + t189 * t192, -t181 * t193, -t174, -t174, 0; t175, -t186 * t212 + t187 * t192, t178 * t193, t172, t172, 0; 0, (t192 * t198 + t193 * t205) * t195, t184 * t193, t176, t176, 0; t178, -t188 * t197, t182, 0, 0, 0; t181, -t186 * t197, -t180, 0, 0, 0; 0, t197 * t209, t185, 0, 0, 0; t172, -t188 * t213 - t189 * t193, -t181 * t192, t175, t175, 0; t174, -t186 * t213 - t187 * t193, t178 * t192, -t216, -t216, 0; 0, (t192 * t205 - t193 * t198) * t195, t184 * t192, t177, t177, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end