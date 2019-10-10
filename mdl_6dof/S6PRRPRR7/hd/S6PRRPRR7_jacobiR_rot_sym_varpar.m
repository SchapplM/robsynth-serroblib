% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:18
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:19
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
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:19
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (19->13), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
	t81 = sin(pkin(6));
	t84 = sin(qJ(3));
	t93 = t81 * t84;
	t85 = sin(qJ(2));
	t92 = t81 * t85;
	t86 = cos(qJ(3));
	t91 = t81 * t86;
	t87 = cos(qJ(2));
	t90 = t81 * t87;
	t83 = cos(pkin(6));
	t89 = t83 * t85;
	t88 = t83 * t87;
	t82 = cos(pkin(11));
	t80 = sin(pkin(11));
	t79 = -t80 * t89 + t82 * t87;
	t78 = -t80 * t88 - t82 * t85;
	t77 = t80 * t87 + t82 * t89;
	t76 = -t80 * t85 + t82 * t88;
	t1 = [0, t79, 0, 0, 0, 0; 0, t77, 0, 0, 0, 0; 0, t92, 0, 0, 0, 0; 0, -t78 * t86, t79 * t84 - t80 * t91, 0, 0, 0; 0, -t76 * t86, t77 * t84 + t82 * t91, 0, 0, 0; 0, -t86 * t90, -t83 * t86 + t84 * t92, 0, 0, 0; 0, t78 * t84, t79 * t86 + t80 * t93, 0, 0, 0; 0, t76 * t84, t77 * t86 - t82 * t93, 0, 0, 0; 0, t84 * t90, t83 * t84 + t85 * t91, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:19
	% EndTime: 2019-10-09 22:37:19
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (51->27), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
	t110 = sin(pkin(6));
	t114 = sin(qJ(3));
	t126 = t110 * t114;
	t117 = cos(qJ(3));
	t125 = t110 * t117;
	t112 = cos(pkin(6));
	t115 = sin(qJ(2));
	t124 = t112 * t115;
	t118 = cos(qJ(2));
	t123 = t112 * t118;
	t113 = sin(qJ(5));
	t122 = t113 * t114;
	t121 = t113 * t118;
	t116 = cos(qJ(5));
	t120 = t114 * t116;
	t119 = t116 * t118;
	t111 = cos(pkin(11));
	t109 = sin(pkin(11));
	t107 = t112 * t114 + t115 * t125;
	t106 = -t112 * t117 + t115 * t126;
	t105 = -t109 * t124 + t111 * t118;
	t104 = t109 * t123 + t111 * t115;
	t103 = t109 * t118 + t111 * t124;
	t102 = t109 * t115 - t111 * t123;
	t101 = t105 * t117 + t109 * t126;
	t100 = t105 * t114 - t109 * t125;
	t99 = t103 * t117 - t111 * t126;
	t98 = t103 * t114 + t111 * t125;
	t1 = [0, -t104 * t122 + t105 * t116, t101 * t113, 0, t100 * t116 - t104 * t113, 0; 0, -t102 * t122 + t103 * t116, t99 * t113, 0, -t102 * t113 + t98 * t116, 0; 0, (t114 * t121 + t115 * t116) * t110, t107 * t113, 0, t106 * t116 + t110 * t121, 0; 0, -t104 * t120 - t105 * t113, t101 * t116, 0, -t100 * t113 - t104 * t116, 0; 0, -t102 * t120 - t103 * t113, t99 * t116, 0, -t102 * t116 - t98 * t113, 0; 0, (-t113 * t115 + t114 * t119) * t110, t107 * t116, 0, -t106 * t113 + t110 * t119, 0; 0, -t104 * t117, -t100, 0, 0, 0; 0, -t102 * t117, -t98, 0, 0, 0; 0, t118 * t125, -t106, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:19
	% EndTime: 2019-10-09 22:37:19
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (113->28), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->36)
	t135 = qJ(5) + qJ(6);
	t133 = sin(t135);
	t140 = sin(qJ(3));
	t151 = t133 * t140;
	t134 = cos(t135);
	t150 = t134 * t140;
	t137 = sin(pkin(6));
	t149 = t137 * t140;
	t142 = cos(qJ(3));
	t148 = t137 * t142;
	t143 = cos(qJ(2));
	t147 = t137 * t143;
	t139 = cos(pkin(6));
	t141 = sin(qJ(2));
	t146 = t139 * t141;
	t145 = t139 * t143;
	t144 = t140 * t143;
	t138 = cos(pkin(11));
	t136 = sin(pkin(11));
	t131 = t139 * t140 + t141 * t148;
	t130 = -t139 * t142 + t141 * t149;
	t129 = -t136 * t146 + t138 * t143;
	t128 = t136 * t145 + t138 * t141;
	t127 = t136 * t143 + t138 * t146;
	t126 = t136 * t141 - t138 * t145;
	t125 = t129 * t142 + t136 * t149;
	t124 = t129 * t140 - t136 * t148;
	t123 = t127 * t142 - t138 * t149;
	t122 = t127 * t140 + t138 * t148;
	t121 = -t130 * t133 + t134 * t147;
	t120 = t130 * t134 + t133 * t147;
	t119 = -t124 * t133 - t128 * t134;
	t118 = t124 * t134 - t128 * t133;
	t117 = -t122 * t133 - t126 * t134;
	t116 = t122 * t134 - t126 * t133;
	t1 = [0, -t128 * t151 + t129 * t134, t125 * t133, 0, t118, t118; 0, -t126 * t151 + t127 * t134, t123 * t133, 0, t116, t116; 0, (t133 * t144 + t134 * t141) * t137, t131 * t133, 0, t120, t120; 0, -t128 * t150 - t129 * t133, t125 * t134, 0, t119, t119; 0, -t126 * t150 - t127 * t133, t123 * t134, 0, t117, t117; 0, (-t133 * t141 + t134 * t144) * t137, t131 * t134, 0, t121, t121; 0, -t128 * t142, -t124, 0, 0, 0; 0, -t126 * t142, -t122, 0, 0, 0; 0, t142 * t147, -t130, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end