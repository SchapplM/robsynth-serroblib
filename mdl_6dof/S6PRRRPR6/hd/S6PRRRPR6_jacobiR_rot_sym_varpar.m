% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:56
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
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
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
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
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
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
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (51->24), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
	t135 = sin(pkin(6));
	t139 = sin(qJ(3));
	t151 = t135 * t139;
	t142 = cos(qJ(3));
	t150 = t135 * t142;
	t143 = cos(qJ(2));
	t149 = t135 * t143;
	t137 = cos(pkin(6));
	t140 = sin(qJ(2));
	t148 = t137 * t140;
	t147 = t137 * t143;
	t138 = sin(qJ(4));
	t146 = t138 * t142;
	t141 = cos(qJ(4));
	t145 = t141 * t142;
	t144 = t142 * t143;
	t136 = cos(pkin(11));
	t134 = sin(pkin(11));
	t132 = t137 * t139 + t140 * t150;
	t131 = t137 * t142 - t140 * t151;
	t130 = -t134 * t148 + t136 * t143;
	t129 = t134 * t147 + t136 * t140;
	t128 = t134 * t143 + t136 * t148;
	t127 = t134 * t140 - t136 * t147;
	t126 = t130 * t142 + t134 * t151;
	t125 = -t130 * t139 + t134 * t150;
	t124 = t128 * t142 - t136 * t151;
	t123 = -t128 * t139 - t136 * t150;
	t1 = [0, -t129 * t145 + t130 * t138, t125 * t141, -t126 * t138 + t129 * t141, 0, 0; 0, -t127 * t145 + t128 * t138, t123 * t141, -t124 * t138 + t127 * t141, 0, 0; 0, (t138 * t140 + t141 * t144) * t135, t131 * t141, -t132 * t138 - t141 * t149, 0, 0; 0, -t129 * t139, t126, 0, 0, 0; 0, -t127 * t139, t124, 0, 0, 0; 0, t139 * t149, t132, 0, 0, 0; 0, -t129 * t146 - t130 * t141, t125 * t138, t126 * t141 + t129 * t138, 0, 0; 0, -t127 * t146 - t128 * t141, t123 * t138, t124 * t141 + t127 * t138, 0, 0; 0, (t138 * t144 - t140 * t141) * t135, t131 * t138, t132 * t141 - t138 * t149, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:08
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (162->46), mult. (461->91), div. (0->0), fcn. (650->12), ass. (0->51)
	t152 = sin(pkin(6));
	t157 = sin(qJ(3));
	t178 = t152 * t157;
	t161 = cos(qJ(3));
	t177 = t152 * t161;
	t162 = cos(qJ(2));
	t176 = t152 * t162;
	t154 = cos(pkin(6));
	t158 = sin(qJ(2));
	t175 = t154 * t158;
	t174 = t154 * t162;
	t156 = sin(qJ(4));
	t173 = t156 * t161;
	t160 = cos(qJ(4));
	t172 = t160 * t161;
	t171 = t161 * t162;
	t151 = sin(pkin(11));
	t153 = cos(pkin(11));
	t145 = t151 * t162 + t153 * t175;
	t135 = t145 * t161 - t153 * t178;
	t144 = t151 * t158 - t153 * t174;
	t126 = t135 * t156 - t144 * t160;
	t127 = t135 * t160 + t144 * t156;
	t155 = sin(qJ(6));
	t159 = cos(qJ(6));
	t170 = t126 * t159 - t127 * t155;
	t169 = t126 * t155 + t127 * t159;
	t147 = -t151 * t175 + t153 * t162;
	t137 = t147 * t161 + t151 * t178;
	t146 = t151 * t174 + t153 * t158;
	t128 = t137 * t156 - t146 * t160;
	t129 = t137 * t160 + t146 * t156;
	t168 = t128 * t159 - t129 * t155;
	t167 = t128 * t155 + t129 * t159;
	t149 = t154 * t157 + t158 * t177;
	t138 = t149 * t156 + t160 * t176;
	t139 = t149 * t160 - t156 * t176;
	t166 = t138 * t159 - t139 * t155;
	t165 = t138 * t155 + t139 * t159;
	t164 = -t155 * t160 + t156 * t159;
	t163 = t155 * t156 + t159 * t160;
	t148 = t154 * t161 - t158 * t178;
	t141 = (t156 * t158 + t160 * t171) * t152;
	t140 = (t156 * t171 - t158 * t160) * t152;
	t136 = -t147 * t157 + t151 * t177;
	t134 = -t145 * t157 - t153 * t177;
	t133 = -t146 * t172 + t147 * t156;
	t132 = -t146 * t173 - t147 * t160;
	t131 = -t144 * t172 + t145 * t156;
	t130 = -t144 * t173 - t145 * t160;
	t1 = [0, t132 * t155 + t133 * t159, t163 * t136, -t168, 0, t168; 0, t130 * t155 + t131 * t159, t163 * t134, -t170, 0, t170; 0, t140 * t155 + t141 * t159, t163 * t148, -t166, 0, t166; 0, t132 * t159 - t133 * t155, t164 * t136, t167, 0, -t167; 0, t130 * t159 - t131 * t155, t164 * t134, t169, 0, -t169; 0, t140 * t159 - t141 * t155, t164 * t148, t165, 0, -t165; 0, t146 * t157, -t137, 0, 0, 0; 0, t144 * t157, -t135, 0, 0, 0; 0, -t157 * t176, -t149, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end