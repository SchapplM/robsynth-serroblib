% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:38
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPPPR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:58
	% EndTime: 2019-10-24 10:38:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:58
	% EndTime: 2019-10-24 10:38:58
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
	% StartTime: 2019-10-24 10:38:57
	% EndTime: 2019-10-24 10:38:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t6 = qJ(1) + pkin(7);
	t5 = cos(t6);
	t4 = sin(t6);
	t1 = [0, 0, 0, 0, 0; -t5, 0, 0, 0, 0; -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; t4, 0, 0, 0, 0; -t5, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:57
	% EndTime: 2019-10-24 10:38:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->4), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->6)
	t18 = cos(pkin(8));
	t17 = sin(pkin(8));
	t16 = qJ(1) + pkin(7);
	t15 = cos(t16);
	t14 = sin(t16);
	t1 = [0, 0, 0, 0, 0; -t15 * t18, 0, 0, 0, 0; -t14 * t18, 0, 0, 0, 0; 0, 0, 0, 0, 0; t15 * t17, 0, 0, 0, 0; t14 * t17, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t14, 0, 0, 0, 0; t15, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:58
	% EndTime: 2019-10-24 10:38:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->7), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->10)
	t27 = sin(pkin(9));
	t30 = cos(pkin(8));
	t32 = t27 * t30;
	t29 = cos(pkin(9));
	t31 = t29 * t30;
	t28 = sin(pkin(8));
	t26 = qJ(1) + pkin(7);
	t25 = cos(t26);
	t24 = sin(t26);
	t1 = [0, 0, 0, 0, 0; -t24 * t27 - t25 * t31, 0, 0, 0, 0; -t24 * t31 + t25 * t27, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t24 * t29 + t25 * t32, 0, 0, 0, 0; t24 * t32 + t25 * t29, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t25 * t28, 0, 0, 0, 0; -t24 * t28, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:58
	% EndTime: 2019-10-24 10:38:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (48->12), mult. (28->14), div. (0->0), fcn. (48->6), ass. (0->15)
	t51 = qJ(1) + pkin(7);
	t47 = sin(t51);
	t53 = cos(pkin(8));
	t55 = t47 * t53;
	t49 = cos(t51);
	t54 = t49 * t53;
	t52 = sin(pkin(8));
	t50 = pkin(9) + qJ(5);
	t48 = cos(t50);
	t46 = sin(t50);
	t45 = -t47 * t46 - t48 * t54;
	t44 = t46 * t54 - t47 * t48;
	t43 = -t49 * t46 + t48 * t55;
	t42 = t46 * t55 + t49 * t48;
	t1 = [0, 0, 0, 0, -t52 * t46; t45, 0, 0, 0, t42; -t43, 0, 0, 0, -t44; 0, 0, 0, 0, -t52 * t48; t44, 0, 0, 0, t43; t42, 0, 0, 0, t45; 0, 0, 0, 0, 0; -t49 * t52, 0, 0, 0, 0; -t47 * t52, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end