% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRPR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRPR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:20:25
	% EndTime: 2019-12-05 16:20:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:20:25
	% EndTime: 2019-12-05 16:20:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:20:25
	% EndTime: 2019-12-05 16:20:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = pkin(8) + qJ(2);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [0, -t10, 0, 0, 0; 0, t11, 0, 0, 0; 0, 0, 0, 0, 0; 0, -t11, 0, 0, 0; 0, -t10, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:20:25
	% EndTime: 2019-12-05 16:20:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t13 = pkin(8) + qJ(2);
	t11 = sin(t13);
	t14 = sin(qJ(3));
	t19 = t11 * t14;
	t15 = cos(qJ(3));
	t18 = t11 * t15;
	t12 = cos(t13);
	t17 = t12 * t14;
	t16 = t12 * t15;
	t1 = [0, -t18, -t17, 0, 0; 0, t16, -t19, 0, 0; 0, 0, t15, 0, 0; 0, t19, -t16, 0, 0; 0, -t17, -t18, 0, 0; 0, 0, -t14, 0, 0; 0, t12, 0, 0, 0; 0, t11, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:20:25
	% EndTime: 2019-12-05 16:20:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (27->9), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->11)
	t21 = pkin(8) + qJ(2);
	t17 = sin(t21);
	t22 = qJ(3) + pkin(9);
	t18 = sin(t22);
	t26 = t17 * t18;
	t20 = cos(t22);
	t25 = t17 * t20;
	t19 = cos(t21);
	t24 = t19 * t18;
	t23 = t19 * t20;
	t1 = [0, -t25, -t24, 0, 0; 0, t23, -t26, 0, 0; 0, 0, t20, 0, 0; 0, t26, -t23, 0, 0; 0, -t24, -t25, 0, 0; 0, 0, -t18, 0, 0; 0, t19, 0, 0, 0; 0, t17, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:20:25
	% EndTime: 2019-12-05 16:20:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (58->14), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t29 = qJ(3) + pkin(9) + qJ(5);
	t25 = sin(t29);
	t30 = pkin(8) + qJ(2);
	t27 = sin(t30);
	t34 = t27 * t25;
	t26 = cos(t29);
	t33 = t27 * t26;
	t28 = cos(t30);
	t32 = t28 * t25;
	t31 = t28 * t26;
	t1 = [0, -t33, -t32, 0, -t32; 0, t31, -t34, 0, -t34; 0, 0, t26, 0, t26; 0, t34, -t31, 0, -t31; 0, -t32, -t33, 0, -t33; 0, 0, -t25, 0, -t25; 0, t28, 0, 0, 0; 0, t27, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end