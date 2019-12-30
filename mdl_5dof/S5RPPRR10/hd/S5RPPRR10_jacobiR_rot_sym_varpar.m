% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRR10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPPRR10_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR10_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:24:03
	% EndTime: 2019-12-29 16:24:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:24:03
	% EndTime: 2019-12-29 16:24:03
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
	% StartTime: 2019-12-29 16:24:09
	% EndTime: 2019-12-29 16:24:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t10 = cos(qJ(1));
	t9 = sin(qJ(1));
	t8 = cos(pkin(8));
	t7 = sin(pkin(8));
	t1 = [-t9 * t8, 0, 0, 0, 0; t10 * t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9 * t7, 0, 0, 0, 0; -t10 * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0; t10, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:24:03
	% EndTime: 2019-12-29 16:24:03
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t34 = cos(qJ(1));
	t33 = sin(qJ(1));
	t32 = cos(pkin(8));
	t31 = sin(pkin(8));
	t1 = [-t33 * t32, 0, 0, 0, 0; t34 * t32, 0, 0, 0, 0; 0, 0, 0, 0, 0; t34, 0, 0, 0, 0; t33, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t33 * t31, 0, 0, 0, 0; t34 * t31, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:24:03
	% EndTime: 2019-12-29 16:24:03
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (12->10), mult. (36->8), div. (0->0), fcn. (58->6), ass. (0->13)
	t21 = sin(pkin(8));
	t22 = cos(pkin(8));
	t23 = sin(qJ(4));
	t25 = cos(qJ(4));
	t28 = t21 * t25 - t22 * t23;
	t27 = t21 * t23 + t22 * t25;
	t26 = cos(qJ(1));
	t24 = sin(qJ(1));
	t20 = t27 * t26;
	t19 = t28 * t26;
	t18 = t27 * t24;
	t17 = t28 * t24;
	t1 = [-t18, 0, 0, t19, 0; t20, 0, 0, t17, 0; 0, 0, 0, -t27, 0; -t17, 0, 0, -t20, 0; t19, 0, 0, -t18, 0; 0, 0, 0, -t28, 0; -t26, 0, 0, 0, 0; -t24, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:24:03
	% EndTime: 2019-12-29 16:24:03
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (50->17), mult. (56->8), div. (0->0), fcn. (90->6), ass. (0->14)
	t43 = qJ(4) + qJ(5);
	t41 = sin(t43);
	t42 = cos(t43);
	t44 = sin(pkin(8));
	t45 = cos(pkin(8));
	t39 = t45 * t41 - t44 * t42;
	t48 = t44 * t41 + t45 * t42;
	t47 = cos(qJ(1));
	t46 = sin(qJ(1));
	t37 = t48 * t47;
	t36 = t39 * t47;
	t35 = t48 * t46;
	t34 = t39 * t46;
	t1 = [-t35, 0, 0, -t36, -t36; t37, 0, 0, -t34, -t34; 0, 0, 0, -t48, -t48; t34, 0, 0, -t37, -t37; -t36, 0, 0, -t35, -t35; 0, 0, 0, t39, t39; -t47, 0, 0, 0, 0; -t46, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end