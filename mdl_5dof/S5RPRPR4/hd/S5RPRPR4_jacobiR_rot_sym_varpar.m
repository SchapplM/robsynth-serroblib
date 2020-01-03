% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRPR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:40:11
	% EndTime: 2020-01-03 11:40:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:40:11
	% EndTime: 2020-01-03 11:40:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t4, 0, 0, 0, 0; t3, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t3, 0, 0, 0, 0; t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:40:11
	% EndTime: 2020-01-03 11:40:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t6 = qJ(1) + pkin(8);
	t5 = cos(t6);
	t4 = sin(t6);
	t1 = [0, 0, 0, 0, 0; t5, 0, 0, 0, 0; t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t4, 0, 0, 0, 0; t5, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:40:11
	% EndTime: 2020-01-03 11:40:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t25 = qJ(1) + pkin(8);
	t23 = sin(t25);
	t26 = sin(qJ(3));
	t30 = t23 * t26;
	t27 = cos(qJ(3));
	t29 = t23 * t27;
	t24 = cos(t25);
	t28 = t24 * t26;
	t22 = t24 * t27;
	t1 = [0, 0, t27, 0, 0; t22, 0, -t30, 0, 0; t29, 0, t28, 0, 0; 0, 0, -t26, 0, 0; -t28, 0, -t29, 0, 0; -t30, 0, t22, 0, 0; 0, 0, 0, 0, 0; t23, 0, 0, 0, 0; -t24, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:40:11
	% EndTime: 2020-01-03 11:40:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (26->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->11)
	t34 = qJ(3) + pkin(9);
	t30 = sin(t34);
	t35 = qJ(1) + pkin(8);
	t31 = sin(t35);
	t38 = t31 * t30;
	t32 = cos(t34);
	t37 = t31 * t32;
	t33 = cos(t35);
	t36 = t33 * t30;
	t29 = t33 * t32;
	t1 = [0, 0, t32, 0, 0; t29, 0, -t38, 0, 0; t37, 0, t36, 0, 0; 0, 0, -t30, 0, 0; -t36, 0, -t37, 0, 0; -t38, 0, t29, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t33, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:40:11
	% EndTime: 2020-01-03 11:40:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (55->11), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t49 = qJ(3) + pkin(9) + qJ(5);
	t45 = sin(t49);
	t50 = qJ(1) + pkin(8);
	t47 = sin(t50);
	t52 = t47 * t45;
	t46 = cos(t49);
	t51 = t47 * t46;
	t48 = cos(t50);
	t44 = t48 * t46;
	t43 = t48 * t45;
	t1 = [0, 0, t46, 0, t46; t44, 0, -t52, 0, -t52; t51, 0, t43, 0, t43; 0, 0, -t45, 0, -t45; -t43, 0, -t51, 0, -t51; -t52, 0, t44, 0, t44; 0, 0, 0, 0, 0; t47, 0, 0, 0, 0; -t48, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end