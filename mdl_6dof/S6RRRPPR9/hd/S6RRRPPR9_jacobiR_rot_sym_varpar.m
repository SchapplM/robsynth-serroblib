% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:31
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR9_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR9_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
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
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
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
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:25
	% DurationCPUTime: 0.07s
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
	% StartTime: 2019-10-10 11:31:25
	% EndTime: 2019-10-10 11:31:25
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (54->25), mult. (163->57), div. (0->0), fcn. (238->10), ass. (0->30)
	t111 = sin(pkin(11));
	t118 = cos(qJ(3));
	t130 = t111 * t118;
	t112 = sin(pkin(6));
	t115 = sin(qJ(3));
	t129 = t112 * t115;
	t128 = t112 * t118;
	t120 = cos(qJ(1));
	t127 = t112 * t120;
	t113 = cos(pkin(11));
	t126 = t113 * t118;
	t116 = sin(qJ(2));
	t117 = sin(qJ(1));
	t125 = t117 * t116;
	t119 = cos(qJ(2));
	t124 = t117 * t119;
	t123 = t118 * t119;
	t122 = t120 * t116;
	t121 = t120 * t119;
	t114 = cos(pkin(6));
	t107 = t114 * t122 + t124;
	t101 = -t107 * t115 - t118 * t127;
	t102 = -t107 * t118 + t115 * t127;
	t109 = -t114 * t125 + t121;
	t108 = t114 * t124 + t122;
	t106 = t114 * t121 - t125;
	t105 = t114 * t118 - t116 * t129;
	t104 = t109 * t118 + t117 * t129;
	t103 = t109 * t115 - t117 * t128;
	t1 = [t102 * t113 + t106 * t111, -t108 * t126 + t109 * t111, -t103 * t113, 0, 0, 0; t104 * t113 + t108 * t111, t106 * t126 + t107 * t111, t101 * t113, 0, 0, 0; 0, (t111 * t116 + t113 * t123) * t112, t105 * t113, 0, 0, 0; -t102 * t111 + t106 * t113, t108 * t130 + t109 * t113, t103 * t111, 0, 0, 0; -t104 * t111 + t108 * t113, -t106 * t130 + t107 * t113, -t101 * t111, 0, 0, 0; 0, (-t111 * t123 + t113 * t116) * t112, -t105 * t111, 0, 0, 0; t101, -t108 * t115, t104, 0, 0, 0; t103, t106 * t115, -t102, 0, 0, 0; 0, t119 * t129, t114 * t115 + t116 * t128, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:25
	% EndTime: 2019-10-10 11:31:25
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (51->24), mult. (163->57), div. (0->0), fcn. (238->10), ass. (0->30)
	t136 = sin(pkin(11));
	t143 = cos(qJ(3));
	t155 = t136 * t143;
	t137 = sin(pkin(6));
	t140 = sin(qJ(3));
	t154 = t137 * t140;
	t153 = t137 * t143;
	t145 = cos(qJ(1));
	t152 = t137 * t145;
	t138 = cos(pkin(11));
	t151 = t138 * t143;
	t141 = sin(qJ(2));
	t142 = sin(qJ(1));
	t150 = t142 * t141;
	t144 = cos(qJ(2));
	t149 = t142 * t144;
	t148 = t143 * t144;
	t147 = t145 * t141;
	t146 = t145 * t144;
	t139 = cos(pkin(6));
	t132 = t139 * t147 + t149;
	t126 = -t132 * t140 - t143 * t152;
	t127 = -t132 * t143 + t140 * t152;
	t134 = -t139 * t150 + t146;
	t133 = t139 * t149 + t147;
	t131 = t139 * t146 - t150;
	t130 = t139 * t143 - t141 * t154;
	t129 = t134 * t143 + t142 * t154;
	t128 = t134 * t140 - t142 * t153;
	t1 = [t127 * t138 + t131 * t136, -t133 * t151 + t134 * t136, -t128 * t138, 0, 0, 0; t129 * t138 + t133 * t136, t131 * t151 + t132 * t136, t126 * t138, 0, 0, 0; 0, (t136 * t141 + t138 * t148) * t137, t130 * t138, 0, 0, 0; t126, -t133 * t140, t129, 0, 0, 0; t128, t131 * t140, -t127, 0, 0, 0; 0, t144 * t154, t139 * t140 + t141 * t153, 0, 0, 0; t127 * t136 - t131 * t138, -t133 * t155 - t134 * t138, -t128 * t136, 0, 0, 0; t129 * t136 - t133 * t138, t131 * t155 - t132 * t138, t126 * t136, 0, 0, 0; 0, (t136 * t148 - t138 * t141) * t137, t130 * t136, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:25
	% EndTime: 2019-10-10 11:31:25
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (156->44), mult. (445->90), div. (0->0), fcn. (628->12), ass. (0->52)
	t176 = cos(pkin(6));
	t179 = sin(qJ(2));
	t184 = cos(qJ(1));
	t189 = t184 * t179;
	t180 = sin(qJ(1));
	t183 = cos(qJ(2));
	t191 = t180 * t183;
	t168 = t176 * t189 + t191;
	t178 = sin(qJ(3));
	t182 = cos(qJ(3));
	t174 = sin(pkin(6));
	t194 = t174 * t184;
	t159 = t168 * t182 - t178 * t194;
	t188 = t184 * t183;
	t192 = t180 * t179;
	t167 = -t176 * t188 + t192;
	t173 = sin(pkin(11));
	t175 = cos(pkin(11));
	t146 = t159 * t173 - t167 * t175;
	t147 = t159 * t175 + t167 * t173;
	t177 = sin(qJ(6));
	t181 = cos(qJ(6));
	t202 = t146 * t181 - t147 * t177;
	t201 = -t146 * t177 - t147 * t181;
	t198 = t173 * t182;
	t197 = t174 * t178;
	t196 = t174 * t182;
	t195 = t174 * t183;
	t193 = t175 * t182;
	t190 = t182 * t183;
	t187 = t173 * t181 - t175 * t177;
	t186 = t173 * t177 + t175 * t181;
	t185 = t168 * t178 + t182 * t194;
	t170 = -t176 * t192 + t188;
	t169 = t176 * t191 + t189;
	t166 = t176 * t178 + t179 * t196;
	t165 = t176 * t182 - t179 * t197;
	t164 = (t173 * t179 + t175 * t190) * t174;
	t163 = (t173 * t190 - t175 * t179) * t174;
	t162 = t170 * t182 + t180 * t197;
	t161 = -t170 * t178 + t180 * t196;
	t157 = t166 * t175 - t173 * t195;
	t156 = t166 * t173 + t175 * t195;
	t155 = -t169 * t193 + t170 * t173;
	t154 = -t169 * t198 - t170 * t175;
	t153 = -t167 * t193 + t168 * t173;
	t152 = -t167 * t198 - t168 * t175;
	t151 = t162 * t175 + t169 * t173;
	t150 = t162 * t173 - t169 * t175;
	t145 = t150 * t177 + t151 * t181;
	t144 = t150 * t181 - t151 * t177;
	t1 = [t201, t154 * t177 + t155 * t181, t186 * t161, 0, 0, t144; t145, t152 * t177 + t153 * t181, -t186 * t185, 0, 0, t202; 0, t163 * t177 + t164 * t181, t186 * t165, 0, 0, t156 * t181 - t157 * t177; -t202, t154 * t181 - t155 * t177, t187 * t161, 0, 0, -t145; t144, t152 * t181 - t153 * t177, -t187 * t185, 0, 0, t201; 0, t163 * t181 - t164 * t177, t187 * t165, 0, 0, -t156 * t177 - t157 * t181; t185, t169 * t178, -t162, 0, 0, 0; t161, t167 * t178, -t159, 0, 0, 0; 0, -t178 * t195, -t166, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end