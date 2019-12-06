% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPPRR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:42:51
	% EndTime: 2019-12-05 17:42:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:42:51
	% EndTime: 2019-12-05 17:42:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t4, 0, 0, 0, 0; -t3, 0, 0, 0, 0; 0, 0, 0, 0, 0; t3, 0, 0, 0, 0; -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:42:51
	% EndTime: 2019-12-05 17:42:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t6 = qJ(1) + pkin(8);
	t5 = cos(t6);
	t4 = sin(t6);
	t1 = [0, 0, 0, 0, 0; -t5, 0, 0, 0, 0; -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; t4, 0, 0, 0, 0; -t5, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:42:51
	% EndTime: 2019-12-05 17:42:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->4), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->6)
	t18 = cos(pkin(9));
	t17 = sin(pkin(9));
	t16 = qJ(1) + pkin(8);
	t15 = cos(t16);
	t14 = sin(t16);
	t1 = [0, 0, 0, 0, 0; -t15 * t18, 0, 0, 0, 0; -t14 * t18, 0, 0, 0, 0; 0, 0, 0, 0, 0; t15 * t17, 0, 0, 0, 0; t14 * t17, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t14, 0, 0, 0, 0; t15, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:42:51
	% EndTime: 2019-12-05 17:42:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (26->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->11)
	t32 = qJ(1) + pkin(8);
	t28 = sin(t32);
	t31 = pkin(9) + qJ(4);
	t29 = cos(t31);
	t35 = t28 * t29;
	t27 = sin(t31);
	t30 = cos(t32);
	t34 = t30 * t27;
	t33 = t30 * t29;
	t26 = t28 * t27;
	t1 = [0, 0, 0, t29, 0; -t33, 0, 0, t26, 0; -t35, 0, 0, -t34, 0; 0, 0, 0, -t27, 0; t34, 0, 0, t35, 0; t26, 0, 0, -t33, 0; 0, 0, 0, 0, 0; -t28, 0, 0, 0, 0; t30, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:42:51
	% EndTime: 2019-12-05 17:42:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (55->11), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t48 = pkin(9) + qJ(4) + qJ(5);
	t44 = sin(t48);
	t49 = qJ(1) + pkin(8);
	t47 = cos(t49);
	t51 = t47 * t44;
	t45 = cos(t48);
	t50 = t47 * t45;
	t46 = sin(t49);
	t43 = t46 * t45;
	t42 = t46 * t44;
	t1 = [0, 0, 0, t45, t45; -t50, 0, 0, t42, t42; -t43, 0, 0, -t51, -t51; 0, 0, 0, -t44, -t44; t51, 0, 0, t43, t43; t42, 0, 0, -t50, -t50; 0, 0, 0, 0, 0; -t46, 0, 0, 0, 0; t47, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end