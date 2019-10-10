% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
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
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; -t6, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t6, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t11 = sin(qJ(3));
	t12 = sin(qJ(1));
	t16 = t12 * t11;
	t13 = cos(qJ(3));
	t14 = cos(qJ(1));
	t15 = t14 * t13;
	t10 = t14 * t11;
	t9 = t12 * t13;
	t1 = [t10, 0, t9, 0, 0, 0; t16, 0, -t15, 0, 0, 0; 0, 0, -t11, 0, 0, 0; t15, 0, -t16, 0, 0, 0; t9, 0, t10, 0, 0, 0; 0, 0, -t13, 0, 0, 0; -t12, 0, 0, 0, 0, 0; t14, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (25->10), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t24 = qJ(3) + qJ(4);
	t22 = sin(t24);
	t25 = sin(qJ(1));
	t28 = t25 * t22;
	t23 = cos(t24);
	t26 = cos(qJ(1));
	t27 = t26 * t23;
	t21 = t26 * t22;
	t20 = t25 * t23;
	t1 = [t21, 0, t20, t20, 0, 0; t28, 0, -t27, -t27, 0, 0; 0, 0, -t22, -t22, 0, 0; t27, 0, -t28, -t28, 0, 0; t20, 0, t21, t21, 0, 0; 0, 0, -t23, -t23, 0, 0; -t25, 0, 0, 0, 0, 0; t26, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (50->19), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t94 = qJ(3) + qJ(4);
	t92 = sin(t94);
	t97 = cos(qJ(5));
	t106 = t92 * t97;
	t95 = sin(qJ(5));
	t96 = sin(qJ(1));
	t105 = t96 * t95;
	t104 = t96 * t97;
	t98 = cos(qJ(1));
	t103 = t98 * t92;
	t102 = t98 * t95;
	t101 = t98 * t97;
	t93 = cos(t94);
	t100 = t93 * t105;
	t99 = t93 * t101;
	t91 = t96 * t92;
	t90 = t92 * t95;
	t89 = t93 * t102;
	t88 = t93 * t104;
	t87 = t92 * t101 - t105;
	t86 = t92 * t102 + t104;
	t85 = t92 * t104 + t102;
	t84 = -t92 * t105 + t101;
	t1 = [t87, 0, t88, t88, t84, 0; t85, 0, -t99, -t99, t86, 0; 0, 0, -t106, -t106, -t93 * t95, 0; -t86, 0, -t100, -t100, -t85, 0; t84, 0, t89, t89, t87, 0; 0, 0, t90, t90, -t93 * t97, 0; -t98 * t93, 0, t91, t91, 0, 0; -t96 * t93, 0, -t103, -t103, 0, 0; 0, 0, t93, t93, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (102->23), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->27)
	t116 = qJ(3) + qJ(4);
	t112 = sin(t116);
	t115 = qJ(5) + qJ(6);
	t113 = cos(t115);
	t128 = t112 * t113;
	t111 = sin(t115);
	t114 = cos(t116);
	t127 = t114 * t111;
	t126 = t114 * t113;
	t117 = sin(qJ(1));
	t110 = t117 * t112;
	t125 = t117 * t113;
	t124 = t117 * t114;
	t118 = cos(qJ(1));
	t123 = t118 * t112;
	t122 = t118 * t113;
	t121 = t118 * t114;
	t120 = t111 * t124;
	t119 = t113 * t121;
	t109 = t112 * t111;
	t108 = t111 * t121;
	t107 = t113 * t124;
	t106 = -t117 * t111 + t112 * t122;
	t105 = t111 * t123 + t125;
	t104 = t118 * t111 + t112 * t125;
	t103 = -t111 * t110 + t122;
	t1 = [t106, 0, t107, t107, t103, t103; t104, 0, -t119, -t119, t105, t105; 0, 0, -t128, -t128, -t127, -t127; -t105, 0, -t120, -t120, -t104, -t104; t103, 0, t108, t108, t106, t106; 0, 0, t109, t109, -t126, -t126; -t121, 0, t110, t110, 0, 0; -t124, 0, -t123, -t123, 0, 0; 0, 0, t114, t114, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end