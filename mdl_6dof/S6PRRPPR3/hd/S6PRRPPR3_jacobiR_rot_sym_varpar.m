% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPPR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
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
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
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
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (16->10), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
	t84 = sin(pkin(6));
	t87 = sin(qJ(3));
	t96 = t84 * t87;
	t88 = sin(qJ(2));
	t95 = t84 * t88;
	t89 = cos(qJ(3));
	t94 = t84 * t89;
	t90 = cos(qJ(2));
	t93 = t84 * t90;
	t86 = cos(pkin(6));
	t92 = t86 * t88;
	t91 = t86 * t90;
	t85 = cos(pkin(10));
	t83 = sin(pkin(10));
	t82 = -t83 * t92 + t85 * t90;
	t81 = -t83 * t91 - t85 * t88;
	t80 = t83 * t90 + t85 * t92;
	t79 = -t83 * t88 + t85 * t91;
	t1 = [0, t81 * t89, -t82 * t87 + t83 * t94, 0, 0, 0; 0, t79 * t89, -t80 * t87 - t85 * t94, 0, 0, 0; 0, t89 * t93, t86 * t89 - t87 * t95, 0, 0, 0; 0, t82, 0, 0, 0, 0; 0, t80, 0, 0, 0, 0; 0, t95, 0, 0, 0, 0; 0, t81 * t87, t82 * t89 + t83 * t96, 0, 0, 0; 0, t79 * t87, t80 * t89 - t85 * t96, 0, 0, 0; 0, t87 * t93, t86 * t87 + t88 * t94, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (20->16), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
	t76 = sin(pkin(6));
	t79 = sin(qJ(3));
	t88 = t76 * t79;
	t80 = sin(qJ(2));
	t87 = t76 * t80;
	t81 = cos(qJ(3));
	t86 = t76 * t81;
	t82 = cos(qJ(2));
	t85 = t76 * t82;
	t78 = cos(pkin(6));
	t84 = t78 * t80;
	t83 = t78 * t82;
	t77 = cos(pkin(10));
	t75 = sin(pkin(10));
	t74 = -t75 * t84 + t77 * t82;
	t73 = -t75 * t83 - t77 * t80;
	t72 = t75 * t82 + t77 * t84;
	t71 = -t75 * t80 + t77 * t83;
	t1 = [0, t73 * t79, t74 * t81 + t75 * t88, 0, 0, 0; 0, t71 * t79, t72 * t81 - t77 * t88, 0, 0, 0; 0, t79 * t85, t78 * t79 + t80 * t86, 0, 0, 0; 0, -t73 * t81, t74 * t79 - t75 * t86, 0, 0, 0; 0, -t71 * t81, t72 * t79 + t77 * t86, 0, 0, 0; 0, -t81 * t85, -t78 * t81 + t79 * t87, 0, 0, 0; 0, -t74, 0, 0, 0, 0; 0, -t72, 0, 0, 0, 0; 0, -t87, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:55
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (54->28), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
	t109 = sin(pkin(6));
	t113 = sin(qJ(3));
	t125 = t109 * t113;
	t116 = cos(qJ(3));
	t124 = t109 * t116;
	t111 = cos(pkin(6));
	t114 = sin(qJ(2));
	t123 = t111 * t114;
	t117 = cos(qJ(2));
	t122 = t111 * t117;
	t112 = sin(qJ(6));
	t121 = t112 * t113;
	t120 = t112 * t117;
	t115 = cos(qJ(6));
	t119 = t113 * t115;
	t118 = t115 * t117;
	t110 = cos(pkin(10));
	t108 = sin(pkin(10));
	t106 = t111 * t113 + t114 * t124;
	t105 = -t111 * t116 + t114 * t125;
	t104 = -t108 * t123 + t110 * t117;
	t103 = -t108 * t122 - t110 * t114;
	t102 = t108 * t117 + t110 * t123;
	t101 = -t108 * t114 + t110 * t122;
	t100 = t104 * t116 + t108 * t125;
	t99 = t104 * t113 - t108 * t124;
	t98 = t102 * t116 - t110 * t125;
	t97 = t102 * t113 + t110 * t124;
	t1 = [0, t103 * t119 - t104 * t112, t100 * t115, 0, 0, t103 * t115 - t99 * t112; 0, t101 * t119 - t102 * t112, t98 * t115, 0, 0, t101 * t115 - t97 * t112; 0, (-t112 * t114 + t113 * t118) * t109, t106 * t115, 0, 0, -t105 * t112 + t109 * t118; 0, -t103 * t121 - t104 * t115, -t100 * t112, 0, 0, -t103 * t112 - t99 * t115; 0, -t101 * t121 - t102 * t115, -t98 * t112, 0, 0, -t101 * t112 - t97 * t115; 0, (-t113 * t120 - t114 * t115) * t109, -t106 * t112, 0, 0, -t105 * t115 - t109 * t120; 0, t103 * t116, -t99, 0, 0, 0; 0, t101 * t116, -t97, 0, 0, 0; 0, t117 * t124, -t105, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end