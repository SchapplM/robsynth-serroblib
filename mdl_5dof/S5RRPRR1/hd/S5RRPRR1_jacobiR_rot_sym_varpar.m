% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRPRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_jacobiR_rot_sym_varpar: pkin has to be [4x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:26:04
	% EndTime: 2019-12-05 18:26:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:26:04
	% EndTime: 2019-12-05 18:26:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:26:04
	% EndTime: 2019-12-05 18:26:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t7 = sin(qJ(2));
	t8 = sin(qJ(1));
	t14 = t8 * t7;
	t9 = cos(qJ(2));
	t13 = t8 * t9;
	t10 = cos(qJ(1));
	t12 = t10 * t7;
	t11 = t10 * t9;
	t1 = [-t13, -t12, 0, 0, 0; t11, -t14, 0, 0, 0; 0, t9, 0, 0, 0; t14, -t11, 0, 0, 0; -t12, -t13, 0, 0, 0; 0, -t7, 0, 0, 0; t10, 0, 0, 0, 0; t8, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:26:04
	% EndTime: 2019-12-05 18:26:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(2));
	t11 = sin(qJ(1));
	t17 = t11 * t10;
	t12 = cos(qJ(2));
	t16 = t11 * t12;
	t13 = cos(qJ(1));
	t15 = t13 * t10;
	t14 = t13 * t12;
	t1 = [-t16, -t15, 0, 0, 0; t14, -t17, 0, 0, 0; 0, t12, 0, 0, 0; t17, -t14, 0, 0, 0; -t15, -t16, 0, 0, 0; 0, -t10, 0, 0, 0; t13, 0, 0, 0, 0; t11, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:26:04
	% EndTime: 2019-12-05 18:26:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t22 = qJ(2) + qJ(4);
	t20 = sin(t22);
	t23 = sin(qJ(1));
	t28 = t23 * t20;
	t21 = cos(t22);
	t27 = t23 * t21;
	t24 = cos(qJ(1));
	t26 = t24 * t20;
	t25 = t24 * t21;
	t1 = [-t27, -t26, 0, -t26, 0; t25, -t28, 0, -t28, 0; 0, t21, 0, t21, 0; t28, -t25, 0, -t25, 0; -t26, -t27, 0, -t27, 0; 0, -t20, 0, -t20, 0; t24, 0, 0, 0, 0; t23, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:26:05
	% EndTime: 2019-12-05 18:26:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (47->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t100 = sin(qJ(5));
	t99 = qJ(2) + qJ(4);
	t98 = cos(t99);
	t110 = t98 * t100;
	t101 = sin(qJ(1));
	t109 = t101 * t100;
	t102 = cos(qJ(5));
	t108 = t101 * t102;
	t103 = cos(qJ(1));
	t107 = t103 * t100;
	t106 = t103 * t102;
	t97 = sin(t99);
	t105 = t97 * t108;
	t104 = t97 * t106;
	t96 = t103 * t98;
	t95 = t98 * t102;
	t94 = t101 * t98;
	t93 = t97 * t107;
	t92 = t97 * t109;
	t91 = t98 * t106 + t109;
	t90 = -t98 * t107 + t108;
	t89 = -t98 * t108 + t107;
	t88 = t98 * t109 + t106;
	t1 = [t89, -t104, 0, -t104, t90; t91, -t105, 0, -t105, -t88; 0, t95, 0, t95, -t97 * t100; t88, t93, 0, t93, -t91; t90, t92, 0, t92, t89; 0, -t110, 0, -t110, -t97 * t102; -t101 * t97, t96, 0, t96, 0; t103 * t97, t94, 0, t94, 0; 0, t97, 0, t97, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end