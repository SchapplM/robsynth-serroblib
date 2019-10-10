% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:58
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:49
	% EndTime: 2019-10-10 11:58:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:49
	% EndTime: 2019-10-10 11:58:49
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
	% StartTime: 2019-10-10 11:58:49
	% EndTime: 2019-10-10 11:58:49
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t16 = t10 * t9;
	t12 = cos(qJ(1));
	t15 = t12 * t9;
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t1 = [-t14, -t15, 0, 0, 0, 0; t13, -t16, 0, 0, 0, 0; 0, t11, 0, 0, 0, 0; t16, -t13, 0, 0, 0, 0; -t15, -t14, 0, 0, 0, 0; 0, -t9, 0, 0, 0, 0; t12, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:49
	% EndTime: 2019-10-10 11:58:50
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t22 = qJ(2) + qJ(3);
	t20 = sin(t22);
	t23 = sin(qJ(1));
	t28 = t23 * t20;
	t21 = cos(t22);
	t27 = t23 * t21;
	t24 = cos(qJ(1));
	t26 = t24 * t20;
	t25 = t24 * t21;
	t1 = [-t27, -t26, -t26, 0, 0, 0; t25, -t28, -t28, 0, 0, 0; 0, t21, t21, 0, 0, 0; t28, -t25, -t25, 0, 0, 0; -t26, -t27, -t27, 0, 0, 0; 0, -t20, -t20, 0, 0, 0; t24, 0, 0, 0, 0, 0; t23, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:50
	% EndTime: 2019-10-10 11:58:50
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (35->12), mult. (38->18), div. (0->0), fcn. (66->6), ass. (0->20)
	t87 = qJ(2) + qJ(3);
	t86 = cos(t87);
	t88 = sin(pkin(11));
	t98 = t86 * t88;
	t90 = sin(qJ(1));
	t97 = t90 * t88;
	t89 = cos(pkin(11));
	t96 = t90 * t89;
	t91 = cos(qJ(1));
	t95 = t91 * t88;
	t94 = t91 * t89;
	t85 = sin(t87);
	t93 = t85 * t96;
	t92 = t85 * t94;
	t84 = t91 * t86;
	t83 = t90 * t86;
	t82 = t86 * t89;
	t81 = t85 * t95;
	t80 = t85 * t97;
	t1 = [-t86 * t96 + t95, -t92, -t92, 0, 0, 0; t86 * t94 + t97, -t93, -t93, 0, 0, 0; 0, t82, t82, 0, 0, 0; t86 * t97 + t94, t81, t81, 0, 0, 0; -t86 * t95 + t96, t80, t80, 0, 0, 0; 0, -t98, -t98, 0, 0, 0; -t90 * t85, t84, t84, 0, 0, 0; t91 * t85, t83, t83, 0, 0, 0; 0, t85, t85, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:50
	% EndTime: 2019-10-10 11:58:50
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (77->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->23)
	t98 = pkin(11) + qJ(5);
	t94 = sin(t98);
	t99 = qJ(2) + qJ(3);
	t97 = cos(t99);
	t106 = t97 * t94;
	t100 = sin(qJ(1));
	t96 = sin(t99);
	t105 = t100 * t96;
	t92 = t100 * t97;
	t101 = cos(qJ(1));
	t104 = t101 * t96;
	t93 = t101 * t97;
	t95 = cos(t98);
	t103 = t95 * t105;
	t102 = t95 * t104;
	t91 = t97 * t95;
	t90 = t94 * t104;
	t89 = t94 * t105;
	t88 = t100 * t94 + t95 * t93;
	t87 = t100 * t95 - t94 * t93;
	t86 = t101 * t94 - t95 * t92;
	t85 = t101 * t95 + t94 * t92;
	t1 = [t86, -t102, -t102, 0, t87, 0; t88, -t103, -t103, 0, -t85, 0; 0, t91, t91, 0, -t96 * t94, 0; t85, t90, t90, 0, -t88, 0; t87, t89, t89, 0, t86, 0; 0, -t106, -t106, 0, -t96 * t95, 0; -t105, t93, t93, 0, 0, 0; t104, t92, t92, 0, 0, 0; 0, t96, t96, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:50
	% EndTime: 2019-10-10 11:58:50
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (139->21), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->25)
	t115 = pkin(11) + qJ(5) + qJ(6);
	t113 = sin(t115);
	t118 = qJ(2) + qJ(3);
	t116 = sin(t118);
	t127 = t116 * t113;
	t114 = cos(t115);
	t126 = t116 * t114;
	t117 = cos(t118);
	t125 = t117 * t113;
	t119 = sin(qJ(1));
	t124 = t119 * t116;
	t111 = t119 * t117;
	t120 = cos(qJ(1));
	t123 = t120 * t116;
	t112 = t120 * t117;
	t122 = t114 * t124;
	t121 = t114 * t123;
	t110 = t117 * t114;
	t109 = t113 * t123;
	t108 = t113 * t124;
	t107 = t114 * t112 + t119 * t113;
	t106 = -t113 * t112 + t119 * t114;
	t105 = -t114 * t111 + t120 * t113;
	t104 = t113 * t111 + t120 * t114;
	t1 = [t105, -t121, -t121, 0, t106, t106; t107, -t122, -t122, 0, -t104, -t104; 0, t110, t110, 0, -t127, -t127; t104, t109, t109, 0, -t107, -t107; t106, t108, t108, 0, t105, t105; 0, -t125, -t125, 0, -t126, -t126; -t124, t112, t112, 0, 0, 0; t123, t111, t111, 0, 0, 0; 0, t116, t116, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end