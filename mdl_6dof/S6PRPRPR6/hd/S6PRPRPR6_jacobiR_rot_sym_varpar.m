% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t46 = cos(pkin(6));
	t47 = sin(qJ(2));
	t50 = t46 * t47;
	t48 = cos(qJ(2));
	t49 = t46 * t48;
	t45 = cos(pkin(10));
	t44 = sin(pkin(6));
	t43 = sin(pkin(10));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, t43 * t49 + t45 * t47, 0, 0, 0, 0; 0, t43 * t47 - t45 * t49, 0, 0, 0, 0; 0, -t44 * t48, 0, 0, 0, 0; 0, -t43 * t50 + t45 * t48, 0, 0, 0, 0; 0, t43 * t48 + t45 * t50, 0, 0, 0, 0; 0, t44 * t47, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (16->12), mult. (57->31), div. (0->0), fcn. (88->8), ass. (0->18)
	t66 = sin(pkin(6));
	t69 = sin(qJ(4));
	t77 = t66 * t69;
	t71 = cos(qJ(4));
	t76 = t66 * t71;
	t72 = cos(qJ(2));
	t75 = t66 * t72;
	t68 = cos(pkin(6));
	t70 = sin(qJ(2));
	t74 = t68 * t70;
	t73 = t68 * t72;
	t67 = cos(pkin(10));
	t65 = sin(pkin(10));
	t64 = -t65 * t74 + t67 * t72;
	t63 = t65 * t73 + t67 * t70;
	t62 = t65 * t72 + t67 * t74;
	t61 = t65 * t70 - t67 * t73;
	t1 = [0, t64 * t69, 0, t63 * t71 - t65 * t77, 0, 0; 0, t62 * t69, 0, t61 * t71 + t67 * t77, 0, 0; 0, t70 * t77, 0, -t68 * t69 - t71 * t75, 0, 0; 0, t64 * t71, 0, -t63 * t69 - t65 * t76, 0, 0; 0, t62 * t71, 0, -t61 * t69 + t67 * t76, 0, 0; 0, t70 * t76, 0, -t68 * t71 + t69 * t75, 0, 0; 0, -t63, 0, 0, 0, 0; 0, -t61, 0, 0, 0, 0; 0, t75, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (37->22), mult. (107->51), div. (0->0), fcn. (156->10), ass. (0->26)
	t102 = sin(qJ(4));
	t96 = sin(pkin(11));
	t113 = t102 * t96;
	t98 = sin(pkin(6));
	t112 = t102 * t98;
	t99 = cos(pkin(11));
	t111 = t102 * t99;
	t104 = cos(qJ(4));
	t110 = t104 * t98;
	t105 = cos(qJ(2));
	t109 = t105 * t98;
	t101 = cos(pkin(6));
	t103 = sin(qJ(2));
	t108 = t101 * t103;
	t107 = t101 * t105;
	t106 = t102 * t103;
	t100 = cos(pkin(10));
	t97 = sin(pkin(10));
	t94 = -t101 * t102 - t104 * t109;
	t93 = t100 * t105 - t97 * t108;
	t92 = t100 * t103 + t97 * t107;
	t91 = t100 * t108 + t97 * t105;
	t90 = -t100 * t107 + t97 * t103;
	t89 = t100 * t112 + t90 * t104;
	t88 = t92 * t104 - t97 * t112;
	t1 = [0, t93 * t111 - t92 * t96, 0, t88 * t99, 0, 0; 0, t91 * t111 - t90 * t96, 0, t89 * t99, 0, 0; 0, (t105 * t96 + t99 * t106) * t98, 0, t94 * t99, 0, 0; 0, -t93 * t113 - t92 * t99, 0, -t88 * t96, 0, 0; 0, -t91 * t113 - t90 * t99, 0, -t89 * t96, 0, 0; 0, (t105 * t99 - t96 * t106) * t98, 0, -t94 * t96, 0, 0; 0, -t93 * t104, 0, t92 * t102 + t97 * t110, 0, 0; 0, -t91 * t104, 0, -t100 * t110 + t90 * t102, 0, 0; 0, -t103 * t110, 0, t101 * t104 - t102 * t109, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (87->29), mult. (163->64), div. (0->0), fcn. (238->10), ass. (0->31)
	t114 = pkin(11) + qJ(6);
	t112 = sin(t114);
	t119 = sin(qJ(4));
	t131 = t112 * t119;
	t113 = cos(t114);
	t130 = t113 * t119;
	t116 = sin(pkin(6));
	t129 = t116 * t119;
	t120 = sin(qJ(2));
	t128 = t116 * t120;
	t121 = cos(qJ(4));
	t127 = t116 * t121;
	t122 = cos(qJ(2));
	t126 = t116 * t122;
	t118 = cos(pkin(6));
	t125 = t118 * t120;
	t124 = t118 * t122;
	t123 = t119 * t120;
	t117 = cos(pkin(10));
	t115 = sin(pkin(10));
	t110 = t118 * t121 - t119 * t126;
	t109 = -t118 * t119 - t121 * t126;
	t108 = -t115 * t125 + t117 * t122;
	t107 = t115 * t124 + t117 * t120;
	t106 = t115 * t122 + t117 * t125;
	t105 = t115 * t120 - t117 * t124;
	t104 = t105 * t119 - t117 * t127;
	t103 = t105 * t121 + t117 * t129;
	t102 = t107 * t119 + t115 * t127;
	t101 = t107 * t121 - t115 * t129;
	t1 = [0, -t107 * t112 + t108 * t130, 0, t101 * t113, 0, -t102 * t112 + t108 * t113; 0, -t105 * t112 + t106 * t130, 0, t103 * t113, 0, -t104 * t112 + t106 * t113; 0, (t112 * t122 + t113 * t123) * t116, 0, t109 * t113, 0, -t110 * t112 + t113 * t128; 0, -t107 * t113 - t108 * t131, 0, -t101 * t112, 0, -t102 * t113 - t108 * t112; 0, -t105 * t113 - t106 * t131, 0, -t103 * t112, 0, -t104 * t113 - t106 * t112; 0, (-t112 * t123 + t113 * t122) * t116, 0, -t109 * t112, 0, -t110 * t113 - t112 * t128; 0, -t108 * t121, 0, t102, 0, 0; 0, -t106 * t121, 0, t104, 0, 0; 0, -t120 * t127, 0, t110, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end