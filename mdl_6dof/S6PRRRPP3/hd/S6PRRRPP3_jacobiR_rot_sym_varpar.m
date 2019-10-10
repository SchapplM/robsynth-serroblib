% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPP3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPP3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(10));
	t20 = sin(pkin(6));
	t19 = sin(pkin(10));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
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
	t67 = cos(pkin(10));
	t65 = sin(pkin(10));
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
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
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
	t108 = cos(pkin(10));
	t106 = sin(pkin(10));
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
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (54->27), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
	t132 = sin(pkin(6));
	t136 = sin(qJ(3));
	t148 = t132 * t136;
	t139 = cos(qJ(3));
	t147 = t132 * t139;
	t140 = cos(qJ(2));
	t146 = t132 * t140;
	t134 = cos(pkin(6));
	t137 = sin(qJ(2));
	t145 = t134 * t137;
	t144 = t134 * t140;
	t135 = sin(qJ(4));
	t143 = t135 * t139;
	t138 = cos(qJ(4));
	t142 = t138 * t139;
	t141 = t139 * t140;
	t133 = cos(pkin(10));
	t131 = sin(pkin(10));
	t129 = t134 * t136 + t137 * t147;
	t128 = t134 * t139 - t137 * t148;
	t127 = -t131 * t145 + t133 * t140;
	t126 = t131 * t144 + t133 * t137;
	t125 = t131 * t140 + t133 * t145;
	t124 = t131 * t137 - t133 * t144;
	t123 = t127 * t139 + t131 * t148;
	t122 = -t127 * t136 + t131 * t147;
	t121 = t125 * t139 - t133 * t148;
	t120 = -t125 * t136 - t133 * t147;
	t1 = [0, -t126 * t136, t123, 0, 0, 0; 0, -t124 * t136, t121, 0, 0, 0; 0, t136 * t146, t129, 0, 0, 0; 0, t126 * t142 - t127 * t135, -t122 * t138, t123 * t135 - t126 * t138, 0, 0; 0, t124 * t142 - t125 * t135, -t120 * t138, t121 * t135 - t124 * t138, 0, 0; 0, (-t135 * t137 - t138 * t141) * t132, -t128 * t138, t129 * t135 + t138 * t146, 0, 0; 0, -t126 * t143 - t127 * t138, t122 * t135, t123 * t138 + t126 * t135, 0, 0; 0, -t124 * t143 - t125 * t138, t120 * t135, t121 * t138 + t124 * t135, 0, 0; 0, (t135 * t141 - t137 * t138) * t132, t128 * t135, t129 * t138 - t135 * t146, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (51->24), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
	t133 = sin(pkin(6));
	t137 = sin(qJ(3));
	t149 = t133 * t137;
	t140 = cos(qJ(3));
	t148 = t133 * t140;
	t141 = cos(qJ(2));
	t147 = t133 * t141;
	t135 = cos(pkin(6));
	t138 = sin(qJ(2));
	t146 = t135 * t138;
	t145 = t135 * t141;
	t136 = sin(qJ(4));
	t144 = t136 * t140;
	t139 = cos(qJ(4));
	t143 = t139 * t140;
	t142 = t140 * t141;
	t134 = cos(pkin(10));
	t132 = sin(pkin(10));
	t130 = t135 * t137 + t138 * t148;
	t129 = t135 * t140 - t138 * t149;
	t128 = -t132 * t146 + t134 * t141;
	t127 = t132 * t145 + t134 * t138;
	t126 = t132 * t141 + t134 * t146;
	t125 = t132 * t138 - t134 * t145;
	t124 = t128 * t140 + t132 * t149;
	t123 = -t128 * t137 + t132 * t148;
	t122 = t126 * t140 - t134 * t149;
	t121 = -t126 * t137 - t134 * t148;
	t1 = [0, -t127 * t137, t124, 0, 0, 0; 0, -t125 * t137, t122, 0, 0, 0; 0, t137 * t147, t130, 0, 0, 0; 0, -t127 * t144 - t128 * t139, t123 * t136, t124 * t139 + t127 * t136, 0, 0; 0, -t125 * t144 - t126 * t139, t121 * t136, t122 * t139 + t125 * t136, 0, 0; 0, (t136 * t142 - t138 * t139) * t133, t129 * t136, t130 * t139 - t136 * t147, 0, 0; 0, -t127 * t143 + t128 * t136, t123 * t139, -t124 * t136 + t127 * t139, 0, 0; 0, -t125 * t143 + t126 * t136, t121 * t139, -t122 * t136 + t125 * t139, 0, 0; 0, (t136 * t138 + t139 * t142) * t133, t129 * t139, -t130 * t136 - t139 * t147, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end