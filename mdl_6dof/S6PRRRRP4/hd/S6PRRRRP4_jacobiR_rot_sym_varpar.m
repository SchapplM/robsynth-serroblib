% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRP4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:07:40
	% EndTime: 2019-10-09 23:07:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:07:40
	% EndTime: 2019-10-09 23:07:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:07:40
	% EndTime: 2019-10-09 23:07:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(11));
	t20 = sin(pkin(6));
	t19 = sin(pkin(11));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:07:40
	% EndTime: 2019-10-09 23:07:40
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (19->13), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
	t66 = sin(pkin(6));
	t69 = sin(qJ(3));
	t78 = t66 * t69;
	t70 = sin(qJ(2));
	t77 = t66 * t70;
	t71 = cos(qJ(3));
	t76 = t66 * t71;
	t72 = cos(qJ(2));
	t75 = t66 * t72;
	t68 = cos(pkin(6));
	t74 = t68 * t70;
	t73 = t68 * t72;
	t67 = cos(pkin(11));
	t65 = sin(pkin(11));
	t64 = -t65 * t74 + t67 * t72;
	t63 = -t65 * t73 - t67 * t70;
	t62 = t65 * t72 + t67 * t74;
	t61 = -t65 * t70 + t67 * t73;
	t1 = [0, t63 * t71, -t64 * t69 + t65 * t76, 0, 0, 0; 0, t61 * t71, -t62 * t69 - t67 * t76, 0, 0, 0; 0, t71 * t75, t68 * t71 - t69 * t77, 0, 0, 0; 0, -t63 * t69, -t64 * t71 - t65 * t78, 0, 0, 0; 0, -t61 * t69, -t62 * t71 + t67 * t78, 0, 0, 0; 0, -t69 * t75, -t68 * t69 - t70 * t76, 0, 0, 0; 0, t64, 0, 0, 0, 0; 0, t62, 0, 0, 0, 0; 0, t77, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:07:40
	% EndTime: 2019-10-09 23:07:40
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (54->27), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
	t107 = sin(pkin(6));
	t111 = sin(qJ(3));
	t123 = t107 * t111;
	t114 = cos(qJ(3));
	t122 = t107 * t114;
	t115 = cos(qJ(2));
	t121 = t107 * t115;
	t109 = cos(pkin(6));
	t112 = sin(qJ(2));
	t120 = t109 * t112;
	t119 = t109 * t115;
	t110 = sin(qJ(4));
	t118 = t110 * t114;
	t113 = cos(qJ(4));
	t117 = t113 * t114;
	t116 = t114 * t115;
	t108 = cos(pkin(11));
	t106 = sin(pkin(11));
	t104 = t109 * t111 + t112 * t122;
	t103 = t109 * t114 - t112 * t123;
	t102 = -t106 * t120 + t108 * t115;
	t101 = t106 * t119 + t108 * t112;
	t100 = t106 * t115 + t108 * t120;
	t99 = t106 * t112 - t108 * t119;
	t98 = t102 * t114 + t106 * t123;
	t97 = -t102 * t111 + t106 * t122;
	t96 = t100 * t114 - t108 * t123;
	t95 = -t100 * t111 - t108 * t122;
	t1 = [0, -t101 * t117 + t102 * t110, t97 * t113, t101 * t113 - t98 * t110, 0, 0; 0, t100 * t110 - t99 * t117, t95 * t113, -t96 * t110 + t99 * t113, 0, 0; 0, (t110 * t112 + t113 * t116) * t107, t103 * t113, -t104 * t110 - t113 * t121, 0, 0; 0, t101 * t118 + t102 * t113, -t97 * t110, -t101 * t110 - t98 * t113, 0, 0; 0, t100 * t113 + t99 * t118, -t95 * t110, -t99 * t110 - t96 * t113, 0, 0; 0, (-t110 * t116 + t112 * t113) * t107, -t103 * t110, -t104 * t113 + t110 * t121, 0, 0; 0, -t101 * t111, t98, 0, 0, 0; 0, -t99 * t111, t96, 0, 0, 0; 0, t111 * t121, t104, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:07:40
	% EndTime: 2019-10-09 23:07:41
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (116->28), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->36)
	t133 = qJ(4) + qJ(5);
	t131 = sin(t133);
	t140 = cos(qJ(3));
	t149 = t131 * t140;
	t132 = cos(t133);
	t148 = t132 * t140;
	t135 = sin(pkin(6));
	t138 = sin(qJ(3));
	t147 = t135 * t138;
	t146 = t135 * t140;
	t141 = cos(qJ(2));
	t145 = t135 * t141;
	t137 = cos(pkin(6));
	t139 = sin(qJ(2));
	t144 = t137 * t139;
	t143 = t137 * t141;
	t142 = t140 * t141;
	t136 = cos(pkin(11));
	t134 = sin(pkin(11));
	t129 = t137 * t138 + t139 * t146;
	t128 = t137 * t140 - t139 * t147;
	t127 = -t134 * t144 + t136 * t141;
	t126 = t134 * t143 + t136 * t139;
	t125 = t134 * t141 + t136 * t144;
	t124 = t134 * t139 - t136 * t143;
	t123 = t127 * t140 + t134 * t147;
	t122 = -t127 * t138 + t134 * t146;
	t121 = t125 * t140 - t136 * t147;
	t120 = -t125 * t138 - t136 * t146;
	t119 = -t129 * t132 + t131 * t145;
	t118 = -t129 * t131 - t132 * t145;
	t117 = -t123 * t132 - t126 * t131;
	t116 = -t123 * t131 + t126 * t132;
	t115 = -t121 * t132 - t124 * t131;
	t114 = -t121 * t131 + t124 * t132;
	t1 = [0, -t126 * t148 + t127 * t131, t122 * t132, t116, t116, 0; 0, -t124 * t148 + t125 * t131, t120 * t132, t114, t114, 0; 0, (t131 * t139 + t132 * t142) * t135, t128 * t132, t118, t118, 0; 0, t126 * t149 + t127 * t132, -t122 * t131, t117, t117, 0; 0, t124 * t149 + t125 * t132, -t120 * t131, t115, t115, 0; 0, (-t131 * t142 + t132 * t139) * t135, -t128 * t131, t119, t119, 0; 0, -t126 * t138, t123, 0, 0, 0; 0, -t124 * t138, t121, 0, 0, 0; 0, t138 * t145, t129, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:07:41
	% EndTime: 2019-10-09 23:07:41
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (113->25), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->36)
	t168 = qJ(4) + qJ(5);
	t166 = sin(t168);
	t175 = cos(qJ(3));
	t184 = t166 * t175;
	t167 = cos(t168);
	t183 = t167 * t175;
	t170 = sin(pkin(6));
	t173 = sin(qJ(3));
	t182 = t170 * t173;
	t181 = t170 * t175;
	t176 = cos(qJ(2));
	t180 = t170 * t176;
	t172 = cos(pkin(6));
	t174 = sin(qJ(2));
	t179 = t172 * t174;
	t178 = t172 * t176;
	t177 = t175 * t176;
	t171 = cos(pkin(11));
	t169 = sin(pkin(11));
	t164 = t172 * t173 + t174 * t181;
	t163 = t172 * t175 - t174 * t182;
	t162 = -t169 * t179 + t171 * t176;
	t161 = t169 * t178 + t171 * t174;
	t160 = t169 * t176 + t171 * t179;
	t159 = t169 * t174 - t171 * t178;
	t158 = t162 * t175 + t169 * t182;
	t157 = -t162 * t173 + t169 * t181;
	t156 = t160 * t175 - t171 * t182;
	t155 = -t160 * t173 - t171 * t181;
	t154 = t164 * t167 - t166 * t180;
	t153 = -t164 * t166 - t167 * t180;
	t152 = t158 * t167 + t161 * t166;
	t151 = -t158 * t166 + t161 * t167;
	t150 = t156 * t167 + t159 * t166;
	t149 = -t156 * t166 + t159 * t167;
	t1 = [0, -t161 * t183 + t162 * t166, t157 * t167, t151, t151, 0; 0, -t159 * t183 + t160 * t166, t155 * t167, t149, t149, 0; 0, (t166 * t174 + t167 * t177) * t170, t163 * t167, t153, t153, 0; 0, -t161 * t173, t158, 0, 0, 0; 0, -t159 * t173, t156, 0, 0, 0; 0, t173 * t180, t164, 0, 0, 0; 0, -t161 * t184 - t162 * t167, t157 * t166, t152, t152, 0; 0, -t159 * t184 - t160 * t167, t155 * t166, t150, t150, 0; 0, (t166 * t177 - t167 * t174) * t170, t163 * t166, t154, t154, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end