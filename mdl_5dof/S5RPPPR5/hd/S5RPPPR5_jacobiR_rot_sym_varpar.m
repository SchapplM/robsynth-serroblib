% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:51
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPPPR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:51:26
	% EndTime: 2019-12-29 15:51:26
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:51:26
	% EndTime: 2019-12-29 15:51:26
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
	% StartTime: 2019-12-29 15:51:26
	% EndTime: 2019-12-29 15:51:26
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [-t5, 0, 0, 0, 0; t6, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t6, 0, 0, 0, 0; t5, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:51:31
	% EndTime: 2019-12-29 15:51:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->3), mult. (8->4), div. (0->0), fcn. (16->4), ass. (0->7)
	t26 = cos(qJ(1));
	t25 = sin(qJ(1));
	t24 = cos(pkin(7));
	t23 = sin(pkin(7));
	t22 = t25 * t23 + t26 * t24;
	t21 = t26 * t23 - t25 * t24;
	t1 = [t21, 0, 0, 0, 0; t22, 0, 0, 0, 0; 0, 0, 0, 0, 0; t22, 0, 0, 0, 0; -t21, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:51:26
	% EndTime: 2019-12-29 15:51:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (8->4), mult. (16->8), div. (0->0), fcn. (28->6), ass. (0->9)
	t23 = cos(qJ(1));
	t22 = sin(qJ(1));
	t21 = cos(pkin(7));
	t20 = sin(pkin(7));
	t19 = cos(pkin(8));
	t18 = sin(pkin(8));
	t15 = t23 * t20 - t22 * t21;
	t14 = -t22 * t20 - t23 * t21;
	t1 = [t15 * t19, 0, 0, 0, 0; -t14 * t19, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t15 * t18, 0, 0, 0, 0; t14 * t18, 0, 0, 0, 0; 0, 0, 0, 0, 0; t14, 0, 0, 0, 0; t15, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:51:31
	% EndTime: 2019-12-29 15:51:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (27->6), mult. (28->8), div. (0->0), fcn. (50->6), ass. (0->14)
	t39 = cos(qJ(1));
	t38 = sin(qJ(1));
	t32 = sin(pkin(7));
	t33 = cos(pkin(7));
	t24 = -t38 * t32 - t39 * t33;
	t31 = pkin(8) + qJ(5);
	t29 = sin(t31);
	t37 = t24 * t29;
	t30 = cos(t31);
	t36 = t24 * t30;
	t25 = t39 * t32 - t38 * t33;
	t35 = t25 * t29;
	t34 = t25 * t30;
	t1 = [t34, 0, 0, 0, t37; -t36, 0, 0, 0, t35; 0, 0, 0, 0, -t30; -t35, 0, 0, 0, t36; t37, 0, 0, 0, t34; 0, 0, 0, 0, t29; t24, 0, 0, 0, 0; t25, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end