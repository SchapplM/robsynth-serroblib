% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:40
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
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
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.03s
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
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (22->18), mult. (57->31), div. (0->0), fcn. (88->8), ass. (0->18)
	t82 = sin(pkin(6));
	t85 = sin(qJ(4));
	t93 = t82 * t85;
	t87 = cos(qJ(4));
	t92 = t82 * t87;
	t88 = cos(qJ(2));
	t91 = t82 * t88;
	t84 = cos(pkin(6));
	t86 = sin(qJ(2));
	t90 = t84 * t86;
	t89 = t84 * t88;
	t83 = cos(pkin(10));
	t81 = sin(pkin(10));
	t80 = -t81 * t90 + t83 * t88;
	t79 = t81 * t89 + t83 * t86;
	t78 = t81 * t88 + t83 * t90;
	t77 = t81 * t86 - t83 * t89;
	t1 = [0, -t79, 0, 0, 0, 0; 0, -t77, 0, 0, 0, 0; 0, t91, 0, 0, 0, 0; 0, -t80 * t85, 0, -t79 * t87 + t81 * t93, 0, 0; 0, -t78 * t85, 0, -t77 * t87 - t83 * t93, 0, 0; 0, -t86 * t93, 0, t84 * t85 + t87 * t91, 0, 0; 0, -t80 * t87, 0, t79 * t85 + t81 * t92, 0, 0; 0, -t78 * t87, 0, t77 * t85 - t83 * t92, 0, 0; 0, -t86 * t92, 0, t84 * t87 - t85 * t91, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (51->24), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
	t108 = sin(pkin(6));
	t112 = sin(qJ(4));
	t124 = t108 * t112;
	t113 = sin(qJ(2));
	t123 = t108 * t113;
	t115 = cos(qJ(4));
	t122 = t108 * t115;
	t116 = cos(qJ(2));
	t121 = t108 * t116;
	t110 = cos(pkin(6));
	t120 = t110 * t113;
	t119 = t110 * t116;
	t111 = sin(qJ(6));
	t118 = t111 * t115;
	t114 = cos(qJ(6));
	t117 = t114 * t115;
	t109 = cos(pkin(10));
	t107 = sin(pkin(10));
	t105 = t110 * t115 - t112 * t121;
	t104 = t110 * t112 + t115 * t121;
	t103 = -t107 * t120 + t109 * t116;
	t102 = t107 * t119 + t109 * t113;
	t101 = t107 * t116 + t109 * t120;
	t100 = t107 * t113 - t109 * t119;
	t99 = t100 * t112 - t109 * t122;
	t98 = t100 * t115 + t109 * t124;
	t97 = t102 * t112 + t107 * t122;
	t96 = -t102 * t115 + t107 * t124;
	t1 = [0, -t102 * t114 - t103 * t118, 0, t97 * t111, 0, -t103 * t111 + t96 * t114; 0, -t100 * t114 - t101 * t118, 0, t99 * t111, 0, -t101 * t111 - t98 * t114; 0, (-t113 * t118 + t114 * t116) * t108, 0, t105 * t111, 0, t104 * t114 - t111 * t123; 0, t102 * t111 - t103 * t117, 0, t97 * t114, 0, -t103 * t114 - t96 * t111; 0, t100 * t111 - t101 * t117, 0, t99 * t114, 0, -t101 * t114 + t98 * t111; 0, (-t111 * t116 - t113 * t117) * t108, 0, t105 * t114, 0, -t104 * t111 - t114 * t123; 0, t103 * t112, 0, -t96, 0, 0; 0, t101 * t112, 0, t98, 0, 0; 0, t112 * t123, 0, -t104, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end