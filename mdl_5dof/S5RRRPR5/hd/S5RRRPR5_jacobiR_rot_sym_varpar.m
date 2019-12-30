% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRPR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:59:38
	% EndTime: 2019-12-29 19:59:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:59:38
	% EndTime: 2019-12-29 19:59:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:59:38
	% EndTime: 2019-12-29 19:59:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t16 = t10 * t9;
	t12 = cos(qJ(1));
	t15 = t12 * t9;
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t1 = [-t14, -t15, 0, 0, 0; t13, -t16, 0, 0, 0; 0, t11, 0, 0, 0; t16, -t13, 0, 0, 0; -t15, -t14, 0, 0, 0; 0, -t9, 0, 0, 0; t12, 0, 0, 0, 0; t10, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:59:38
	% EndTime: 2019-12-29 19:59:38
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
	t1 = [-t27, -t26, -t26, 0, 0; t25, -t28, -t28, 0, 0; 0, t21, t21, 0, 0; t28, -t25, -t25, 0, 0; -t26, -t27, -t27, 0, 0; 0, -t20, -t20, 0, 0; t24, 0, 0, 0, 0; t23, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:59:38
	% EndTime: 2019-12-29 19:59:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (44->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t27 = qJ(2) + qJ(3) + pkin(9);
	t25 = sin(t27);
	t28 = sin(qJ(1));
	t33 = t28 * t25;
	t26 = cos(t27);
	t32 = t28 * t26;
	t29 = cos(qJ(1));
	t31 = t29 * t25;
	t30 = t29 * t26;
	t1 = [-t32, -t31, -t31, 0, 0; t30, -t33, -t33, 0, 0; 0, t26, t26, 0, 0; t33, -t30, -t30, 0, 0; -t31, -t32, -t32, 0, 0; 0, -t25, -t25, 0, 0; t29, 0, 0, 0, 0; t28, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:59:38
	% EndTime: 2019-12-29 19:59:39
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (77->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t106 = qJ(2) + qJ(3) + pkin(9);
	t105 = cos(t106);
	t107 = sin(qJ(5));
	t117 = t105 * t107;
	t108 = sin(qJ(1));
	t116 = t108 * t107;
	t109 = cos(qJ(5));
	t115 = t108 * t109;
	t110 = cos(qJ(1));
	t114 = t110 * t107;
	t113 = t110 * t109;
	t104 = sin(t106);
	t112 = t104 * t115;
	t111 = t104 * t113;
	t103 = t110 * t105;
	t102 = t105 * t109;
	t101 = t108 * t105;
	t100 = t104 * t114;
	t99 = t104 * t116;
	t98 = t105 * t113 + t116;
	t97 = -t105 * t114 + t115;
	t96 = -t105 * t115 + t114;
	t95 = t105 * t116 + t113;
	t1 = [t96, -t111, -t111, 0, t97; t98, -t112, -t112, 0, -t95; 0, t102, t102, 0, -t104 * t107; t95, t100, t100, 0, -t98; t97, t99, t99, 0, t96; 0, -t117, -t117, 0, -t104 * t109; -t108 * t104, t103, t103, 0, 0; t110 * t104, t101, t101, 0, 0; 0, t104, t104, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end