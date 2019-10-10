% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR9_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR9_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:09
	% EndTime: 2019-10-10 12:44:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:10
	% EndTime: 2019-10-10 12:44:10
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
	% StartTime: 2019-10-10 12:44:10
	% EndTime: 2019-10-10 12:44:10
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
	% StartTime: 2019-10-10 12:44:10
	% EndTime: 2019-10-10 12:44:10
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
	% StartTime: 2019-10-10 12:44:10
	% EndTime: 2019-10-10 12:44:10
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (77->18), mult. (117->30), div. (0->0), fcn. (180->8), ass. (0->28)
	t108 = sin(pkin(6));
	t110 = sin(qJ(2));
	t122 = t108 * t110;
	t111 = sin(qJ(1));
	t121 = t108 * t111;
	t112 = cos(qJ(2));
	t120 = t108 * t112;
	t113 = cos(qJ(1));
	t119 = t108 * t113;
	t118 = t111 * t110;
	t117 = t111 * t112;
	t116 = t113 * t110;
	t115 = t113 * t112;
	t109 = cos(pkin(6));
	t101 = t109 * t116 + t117;
	t107 = qJ(3) + qJ(4);
	t105 = sin(t107);
	t106 = cos(t107);
	t95 = -t101 * t106 + t105 * t119;
	t114 = t101 * t105 + t106 * t119;
	t103 = -t109 * t118 + t115;
	t102 = t109 * t117 + t116;
	t100 = t109 * t115 - t118;
	t99 = -t105 * t109 - t106 * t122;
	t98 = -t105 * t122 + t106 * t109;
	t97 = t103 * t106 + t105 * t121;
	t96 = -t103 * t105 + t106 * t121;
	t1 = [t95, -t102 * t106, t96, t96, 0, 0; t97, t100 * t106, -t114, -t114, 0, 0; 0, t106 * t120, t98, t98, 0, 0; t114, t102 * t105, -t97, -t97, 0, 0; t96, -t100 * t105, t95, t95, 0, 0; 0, -t105 * t120, t99, t99, 0, 0; t100, t103, 0, 0, 0, 0; t102, t101, 0, 0, 0, 0; 0, t122, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:10
	% EndTime: 2019-10-10 12:44:10
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (129->30), mult. (214->58), div. (0->0), fcn. (313->10), ass. (0->38)
	t149 = cos(pkin(6));
	t150 = sin(qJ(2));
	t153 = cos(qJ(1));
	t155 = t153 * t150;
	t151 = sin(qJ(1));
	t152 = cos(qJ(2));
	t156 = t151 * t152;
	t139 = t149 * t155 + t156;
	t145 = qJ(3) + qJ(4);
	t143 = sin(t145);
	t144 = cos(t145);
	t147 = sin(pkin(6));
	t158 = t147 * t153;
	t130 = -t139 * t143 - t144 * t158;
	t146 = sin(pkin(12));
	t166 = t130 * t146;
	t154 = t153 * t152;
	t157 = t151 * t150;
	t141 = -t149 * t157 + t154;
	t159 = t147 * t151;
	t133 = t141 * t143 - t144 * t159;
	t165 = t133 * t146;
	t160 = t147 * t150;
	t136 = -t143 * t160 + t149 * t144;
	t164 = t136 * t146;
	t163 = t144 * t146;
	t148 = cos(pkin(12));
	t162 = t144 * t148;
	t161 = t144 * t152;
	t132 = -t139 * t144 + t143 * t158;
	t140 = t149 * t156 + t155;
	t138 = t149 * t154 - t157;
	t137 = t149 * t143 + t144 * t160;
	t135 = t136 * t148;
	t134 = t141 * t144 + t143 * t159;
	t129 = t133 * t148;
	t128 = t130 * t148;
	t1 = [t132 * t148 + t138 * t146, -t140 * t162 + t141 * t146, -t129, -t129, 0, 0; t134 * t148 + t140 * t146, t138 * t162 + t139 * t146, t128, t128, 0, 0; 0, (t146 * t150 + t148 * t161) * t147, t135, t135, 0, 0; -t132 * t146 + t138 * t148, t140 * t163 + t141 * t148, t165, t165, 0, 0; -t134 * t146 + t140 * t148, -t138 * t163 + t139 * t148, -t166, -t166, 0, 0; 0, (-t146 * t161 + t148 * t150) * t147, -t164, -t164, 0, 0; t130, -t140 * t143, t134, t134, 0, 0; t133, t138 * t143, -t132, -t132, 0, 0; 0, t147 * t152 * t143, t137, t137, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:10
	% EndTime: 2019-10-10 12:44:10
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (205->36), mult. (270->62), div. (0->0), fcn. (395->10), ass. (0->44)
	t173 = cos(pkin(6));
	t174 = sin(qJ(2));
	t177 = cos(qJ(1));
	t179 = t177 * t174;
	t175 = sin(qJ(1));
	t176 = cos(qJ(2));
	t180 = t175 * t176;
	t161 = t173 * t179 + t180;
	t171 = qJ(3) + qJ(4);
	t168 = sin(t171);
	t169 = cos(t171);
	t172 = sin(pkin(6));
	t182 = t172 * t177;
	t154 = -t161 * t169 + t168 * t182;
	t178 = t177 * t176;
	t181 = t175 * t174;
	t160 = -t173 * t178 + t181;
	t170 = pkin(12) + qJ(6);
	t166 = sin(t170);
	t167 = cos(t170);
	t195 = t154 * t166 + t160 * t167;
	t194 = t154 * t167 - t160 * t166;
	t152 = -t161 * t168 - t169 * t182;
	t193 = t152 * t166;
	t163 = -t173 * t181 + t178;
	t184 = t172 * t175;
	t155 = t163 * t168 - t169 * t184;
	t192 = t155 * t166;
	t185 = t172 * t174;
	t158 = -t168 * t185 + t173 * t169;
	t191 = t158 * t166;
	t188 = t166 * t169;
	t187 = t167 * t169;
	t186 = t169 * t176;
	t183 = t172 * t176;
	t162 = t173 * t180 + t179;
	t159 = t173 * t168 + t169 * t185;
	t157 = t158 * t167;
	t156 = t163 * t169 + t168 * t184;
	t151 = t155 * t167;
	t150 = t152 * t167;
	t149 = t156 * t167 + t162 * t166;
	t148 = -t156 * t166 + t162 * t167;
	t1 = [t194, -t162 * t187 + t163 * t166, -t151, -t151, 0, t148; t149, -t160 * t187 + t161 * t166, t150, t150, 0, t195; 0, (t166 * t174 + t167 * t186) * t172, t157, t157, 0, -t159 * t166 - t167 * t183; -t195, t162 * t188 + t163 * t167, t192, t192, 0, -t149; t148, t160 * t188 + t161 * t167, -t193, -t193, 0, t194; 0, (-t166 * t186 + t167 * t174) * t172, -t191, -t191, 0, -t159 * t167 + t166 * t183; t152, -t162 * t168, t156, t156, 0, 0; t155, -t160 * t168, -t154, -t154, 0, 0; 0, t168 * t183, t159, t159, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end