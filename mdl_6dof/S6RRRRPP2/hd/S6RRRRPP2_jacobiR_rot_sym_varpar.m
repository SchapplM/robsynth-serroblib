% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:22
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPP2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:15
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
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:15
	% DurationCPUTime: 0.03s
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
	% StartTime: 2019-10-10 12:22:15
	% EndTime: 2019-10-10 12:22:15
	% DurationCPUTime: 0.06s
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
	% StartTime: 2019-10-10 12:22:16
	% EndTime: 2019-10-10 12:22:16
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (47->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t101 = sin(qJ(4));
	t100 = qJ(2) + qJ(3);
	t99 = cos(t100);
	t111 = t99 * t101;
	t102 = sin(qJ(1));
	t110 = t102 * t101;
	t103 = cos(qJ(4));
	t109 = t102 * t103;
	t104 = cos(qJ(1));
	t108 = t104 * t101;
	t107 = t104 * t103;
	t98 = sin(t100);
	t106 = t98 * t109;
	t105 = t98 * t107;
	t97 = t104 * t99;
	t96 = t99 * t103;
	t95 = t102 * t99;
	t94 = t98 * t108;
	t93 = t98 * t110;
	t92 = t99 * t107 + t110;
	t91 = -t99 * t108 + t109;
	t90 = -t99 * t109 + t108;
	t89 = t99 * t110 + t107;
	t1 = [t90, -t105, -t105, t91, 0, 0; t92, -t106, -t106, -t89, 0, 0; 0, t96, t96, -t98 * t101, 0, 0; t89, t94, t94, -t92, 0, 0; t91, t93, t93, t90, 0, 0; 0, -t111, -t111, -t98 * t103, 0, 0; -t102 * t98, t97, t97, 0, 0, 0; t104 * t98, t95, t95, 0, 0, 0; 0, t98, t98, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:16
	% EndTime: 2019-10-10 12:22:16
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (48->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t115 = sin(qJ(4));
	t116 = sin(qJ(1));
	t126 = t116 * t115;
	t117 = cos(qJ(4));
	t125 = t116 * t117;
	t118 = cos(qJ(1));
	t124 = t118 * t115;
	t123 = t118 * t117;
	t114 = qJ(2) + qJ(3);
	t112 = sin(t114);
	t122 = t112 * t126;
	t121 = t112 * t125;
	t120 = t112 * t124;
	t119 = t112 * t123;
	t113 = cos(t114);
	t111 = t118 * t113;
	t110 = t113 * t117;
	t109 = t113 * t115;
	t108 = t116 * t113;
	t107 = t113 * t123 + t126;
	t106 = t113 * t124 - t125;
	t105 = t113 * t125 - t124;
	t104 = -t113 * t126 - t123;
	t1 = [-t105, -t119, -t119, -t106, 0, 0; t107, -t121, -t121, t104, 0, 0; 0, t110, t110, -t112 * t115, 0, 0; -t116 * t112, t111, t111, 0, 0, 0; t118 * t112, t108, t108, 0, 0, 0; 0, t112, t112, 0, 0, 0; t104, -t120, -t120, t107, 0, 0; t106, -t122, -t122, t105, 0, 0; 0, t109, t109, t112 * t117, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:22:16
	% EndTime: 2019-10-10 12:22:16
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (54->23), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t100 = sin(qJ(1));
	t98 = qJ(2) + qJ(3);
	t97 = cos(t98);
	t112 = t100 * t97;
	t99 = sin(qJ(4));
	t111 = t100 * t99;
	t102 = cos(qJ(1));
	t110 = t102 * t97;
	t109 = t102 * t99;
	t101 = cos(qJ(4));
	t108 = t100 * t101;
	t107 = t102 * t101;
	t96 = sin(t98);
	t106 = t96 * t111;
	t105 = t96 * t109;
	t104 = t96 * t108;
	t103 = t96 * t107;
	t95 = t97 * t101;
	t94 = t97 * t99;
	t93 = t97 * t107 + t111;
	t92 = t97 * t109 - t108;
	t91 = t97 * t108 - t109;
	t90 = -t97 * t111 - t107;
	t1 = [-t91, -t103, -t103, -t92, 0, 0; t93, -t104, -t104, t90, 0, 0; 0, t95, t95, -t96 * t99, 0, 0; t90, -t105, -t105, t93, 0, 0; t92, -t106, -t106, t91, 0, 0; 0, t94, t94, t96 * t101, 0, 0; t100 * t96, -t110, -t110, 0, 0, 0; -t102 * t96, -t112, -t112, 0, 0, 0; 0, -t96, -t96, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end