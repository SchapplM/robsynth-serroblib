% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:17
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:28
	% EndTime: 2019-10-09 23:17:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:28
	% EndTime: 2019-10-09 23:17:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:28
	% EndTime: 2019-10-09 23:17:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(12));
	t20 = sin(pkin(6));
	t19 = sin(pkin(12));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:29
	% EndTime: 2019-10-09 23:17:29
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
	t67 = cos(pkin(12));
	t65 = sin(pkin(12));
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
	% StartTime: 2019-10-09 23:17:29
	% EndTime: 2019-10-09 23:17:29
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
	t108 = cos(pkin(12));
	t106 = sin(pkin(12));
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
	% StartTime: 2019-10-09 23:17:29
	% EndTime: 2019-10-09 23:17:29
	% DurationCPUTime: 0.18s
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
	t136 = cos(pkin(12));
	t134 = sin(pkin(12));
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
	% StartTime: 2019-10-09 23:17:29
	% EndTime: 2019-10-09 23:17:29
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (202->28), mult. (275->63), div. (0->0), fcn. (402->10), ass. (0->36)
	t138 = qJ(4) + qJ(5) + qJ(6);
	t136 = sin(t138);
	t145 = cos(qJ(3));
	t154 = t136 * t145;
	t137 = cos(t138);
	t153 = t137 * t145;
	t140 = sin(pkin(6));
	t143 = sin(qJ(3));
	t152 = t140 * t143;
	t151 = t140 * t145;
	t146 = cos(qJ(2));
	t150 = t140 * t146;
	t142 = cos(pkin(6));
	t144 = sin(qJ(2));
	t149 = t142 * t144;
	t148 = t142 * t146;
	t147 = t145 * t146;
	t141 = cos(pkin(12));
	t139 = sin(pkin(12));
	t134 = t142 * t143 + t144 * t151;
	t133 = t142 * t145 - t144 * t152;
	t132 = -t139 * t149 + t141 * t146;
	t131 = t139 * t148 + t141 * t144;
	t130 = t139 * t146 + t141 * t149;
	t129 = t139 * t144 - t141 * t148;
	t128 = t132 * t145 + t139 * t152;
	t127 = -t132 * t143 + t139 * t151;
	t126 = t130 * t145 - t141 * t152;
	t125 = -t130 * t143 - t141 * t151;
	t124 = -t134 * t137 + t136 * t150;
	t123 = -t134 * t136 - t137 * t150;
	t122 = -t128 * t137 - t131 * t136;
	t121 = -t128 * t136 + t131 * t137;
	t120 = -t126 * t137 - t129 * t136;
	t119 = -t126 * t136 + t129 * t137;
	t1 = [0, -t131 * t153 + t132 * t136, t127 * t137, t121, t121, t121; 0, -t129 * t153 + t130 * t136, t125 * t137, t119, t119, t119; 0, (t136 * t144 + t137 * t147) * t140, t133 * t137, t123, t123, t123; 0, t131 * t154 + t132 * t137, -t127 * t136, t122, t122, t122; 0, t129 * t154 + t130 * t137, -t125 * t136, t120, t120, t120; 0, (-t136 * t147 + t137 * t144) * t140, -t133 * t136, t124, t124, t124; 0, -t131 * t143, t128, 0, 0, 0; 0, -t129 * t143, t126, 0, 0, 0; 0, t143 * t150, t134, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end