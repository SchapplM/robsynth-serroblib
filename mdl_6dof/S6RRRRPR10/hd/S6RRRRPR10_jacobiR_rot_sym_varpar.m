% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR10_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR10_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
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
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
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
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
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
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:05
	% DurationCPUTime: 0.16s
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
	t99 = -t109 * t105 - t106 * t122;
	t98 = -t105 * t122 + t109 * t106;
	t97 = t103 * t106 + t105 * t121;
	t96 = -t103 * t105 + t106 * t121;
	t1 = [t95, -t102 * t106, t96, t96, 0, 0; t97, t100 * t106, -t114, -t114, 0, 0; 0, t106 * t120, t98, t98, 0, 0; t114, t102 * t105, -t97, -t97, 0, 0; t96, -t100 * t105, t95, t95, 0, 0; 0, -t105 * t120, t99, t99, 0, 0; t100, t103, 0, 0, 0, 0; t102, t101, 0, 0, 0, 0; 0, t122, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:05
	% EndTime: 2019-10-10 12:46:06
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (77->16), mult. (117->30), div. (0->0), fcn. (180->8), ass. (0->28)
	t130 = sin(pkin(6));
	t132 = sin(qJ(2));
	t143 = t130 * t132;
	t133 = sin(qJ(1));
	t142 = t130 * t133;
	t134 = cos(qJ(2));
	t141 = t130 * t134;
	t135 = cos(qJ(1));
	t140 = t130 * t135;
	t139 = t133 * t132;
	t138 = t133 * t134;
	t137 = t135 * t132;
	t136 = t135 * t134;
	t131 = cos(pkin(6));
	t124 = t131 * t137 + t138;
	t129 = qJ(3) + qJ(4);
	t127 = sin(t129);
	t128 = cos(t129);
	t117 = t124 * t127 + t128 * t140;
	t118 = t124 * t128 - t127 * t140;
	t126 = -t131 * t139 + t136;
	t125 = t131 * t138 + t137;
	t123 = t131 * t136 - t139;
	t122 = t131 * t127 + t128 * t143;
	t121 = t127 * t143 - t131 * t128;
	t120 = t126 * t128 + t127 * t142;
	t119 = t126 * t127 - t128 * t142;
	t1 = [t123, t126, 0, 0, 0, 0; t125, t124, 0, 0, 0, 0; 0, t143, 0, 0, 0, 0; t118, t125 * t128, t119, t119, 0, 0; -t120, -t123 * t128, t117, t117, 0, 0; 0, -t128 * t141, t121, t121, 0, 0; -t117, -t125 * t127, t120, t120, 0, 0; t119, t123 * t127, t118, t118, 0, 0; 0, t127 * t141, t122, t122, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:46:06
	% EndTime: 2019-10-10 12:46:06
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (155->37), mult. (270->65), div. (0->0), fcn. (395->10), ass. (0->41)
	t169 = sin(qJ(2));
	t170 = sin(qJ(1));
	t172 = cos(qJ(2));
	t173 = cos(qJ(1));
	t186 = cos(pkin(6));
	t175 = t173 * t186;
	t158 = t169 * t175 + t170 * t172;
	t166 = qJ(3) + qJ(4);
	t164 = sin(t166);
	t165 = cos(t166);
	t167 = sin(pkin(6));
	t179 = t167 * t173;
	t148 = t158 * t164 + t165 * t179;
	t157 = t170 * t169 - t172 * t175;
	t168 = sin(qJ(6));
	t171 = cos(qJ(6));
	t188 = -t148 * t168 - t157 * t171;
	t187 = t148 * t171 - t157 * t168;
	t183 = t164 * t168;
	t182 = t164 * t171;
	t181 = t167 * t169;
	t180 = t167 * t170;
	t178 = t168 * t172;
	t177 = t171 * t172;
	t176 = t170 * t186;
	t174 = -t158 * t165 + t164 * t179;
	t160 = -t169 * t176 + t173 * t172;
	t159 = t173 * t169 + t172 * t176;
	t156 = t186 * t164 + t165 * t181;
	t155 = t164 * t181 - t186 * t165;
	t154 = t156 * t171;
	t153 = t156 * t168;
	t152 = t160 * t165 + t164 * t180;
	t151 = t160 * t164 - t165 * t180;
	t147 = t152 * t171;
	t146 = t152 * t168;
	t145 = t174 * t171;
	t144 = t174 * t168;
	t143 = t151 * t168 + t159 * t171;
	t142 = t151 * t171 - t159 * t168;
	t1 = [t188, -t159 * t183 + t160 * t171, t146, t146, 0, t142; t143, -t157 * t183 + t158 * t171, -t144, -t144, 0, t187; 0, (t164 * t178 + t169 * t171) * t167, t153, t153, 0, t155 * t171 + t167 * t178; -t187, -t159 * t182 - t160 * t168, t147, t147, 0, -t143; t142, -t157 * t182 - t158 * t168, -t145, -t145, 0, t188; 0, (t164 * t177 - t168 * t169) * t167, t154, t154, 0, -t155 * t168 + t167 * t177; t174, -t159 * t165, -t151, -t151, 0, 0; t152, -t157 * t165, -t148, -t148, 0, 0; 0, t167 * t172 * t165, -t155, -t155, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end