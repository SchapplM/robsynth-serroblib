% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:45
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
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
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(10);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [-t10, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t11, 0, 0, 0, 0, 0; -t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t13 = qJ(1) + pkin(10);
	t11 = sin(t13);
	t14 = sin(qJ(3));
	t19 = t11 * t14;
	t15 = cos(qJ(3));
	t18 = t11 * t15;
	t12 = cos(t13);
	t17 = t12 * t14;
	t16 = t12 * t15;
	t1 = [-t18, 0, -t17, 0, 0, 0; t16, 0, -t19, 0, 0, 0; 0, 0, t15, 0, 0, 0; t19, 0, -t16, 0, 0, 0; -t17, 0, -t18, 0, 0, 0; 0, 0, -t14, 0, 0, 0; t12, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (42->14), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t26 = qJ(1) + pkin(10);
	t22 = sin(t26);
	t27 = qJ(3) + qJ(4);
	t24 = sin(t27);
	t31 = t22 * t24;
	t25 = cos(t27);
	t30 = t22 * t25;
	t23 = cos(t26);
	t29 = t23 * t24;
	t28 = t23 * t25;
	t1 = [-t30, 0, -t29, -t29, 0, 0; t28, 0, -t31, -t31, 0, 0; 0, 0, t25, t25, 0, 0; t31, 0, -t28, -t28, 0, 0; -t29, 0, -t30, -t30, 0, 0; 0, 0, -t24, -t24, 0, 0; t23, 0, 0, 0, 0, 0; t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (77->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->23)
	t102 = sin(qJ(5));
	t101 = qJ(3) + qJ(4);
	t98 = sin(t101);
	t108 = t98 * t102;
	t103 = cos(qJ(5));
	t107 = t98 * t103;
	t99 = cos(t101);
	t106 = t99 * t102;
	t95 = t99 * t103;
	t100 = qJ(1) + pkin(10);
	t96 = sin(t100);
	t105 = t96 * t107;
	t97 = cos(t100);
	t104 = t97 * t107;
	t94 = t97 * t99;
	t93 = t96 * t99;
	t92 = t97 * t108;
	t91 = t96 * t108;
	t90 = t96 * t102 + t97 * t95;
	t89 = t96 * t103 - t97 * t106;
	t88 = t97 * t102 - t96 * t95;
	t87 = t97 * t103 + t96 * t106;
	t1 = [t88, 0, -t104, -t104, t89, 0; t90, 0, -t105, -t105, -t87, 0; 0, 0, t95, t95, -t108, 0; t87, 0, t92, t92, -t90, 0; t89, 0, t91, t91, t88, 0; 0, 0, -t106, -t106, -t107, 0; -t96 * t98, 0, t94, t94, 0, 0; t97 * t98, 0, t93, t93, 0, 0; 0, 0, t98, t98, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:45:23
	% EndTime: 2019-10-10 01:45:23
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (78->18), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->23)
	t112 = qJ(3) + qJ(4);
	t109 = sin(t112);
	t113 = sin(qJ(5));
	t120 = t109 * t113;
	t114 = cos(qJ(5));
	t119 = t109 * t114;
	t110 = cos(t112);
	t105 = t110 * t113;
	t106 = t110 * t114;
	t111 = qJ(1) + pkin(10);
	t107 = sin(t111);
	t118 = t107 * t120;
	t117 = t107 * t119;
	t108 = cos(t111);
	t116 = t108 * t120;
	t115 = t108 * t119;
	t104 = t108 * t110;
	t103 = t107 * t110;
	t102 = t108 * t106 + t107 * t113;
	t101 = t108 * t105 - t107 * t114;
	t100 = t107 * t106 - t108 * t113;
	t99 = -t107 * t105 - t108 * t114;
	t1 = [-t100, 0, -t115, -t115, -t101, 0; t102, 0, -t117, -t117, t99, 0; 0, 0, t106, t106, -t120, 0; -t107 * t109, 0, t104, t104, 0, 0; t108 * t109, 0, t103, t103, 0, 0; 0, 0, t109, t109, 0, 0; t99, 0, -t116, -t116, t102, 0; t101, 0, -t118, -t118, t100, 0; 0, 0, t105, t105, t119, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end