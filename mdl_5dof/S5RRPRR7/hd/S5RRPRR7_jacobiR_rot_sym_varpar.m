% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRPRR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:04:53
	% EndTime: 2019-12-29 19:04:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:04:51
	% EndTime: 2019-12-29 19:04:51
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
	% StartTime: 2019-12-29 19:04:53
	% EndTime: 2019-12-29 19:04:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (14->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t16 = qJ(1) + qJ(2);
	t15 = cos(t16);
	t14 = sin(t16);
	t1 = [-t14, -t14, 0, 0, 0; t15, t15, 0, 0, 0; 0, 0, 0, 0, 0; -t15, -t15, 0, 0, 0; -t14, -t14, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:04:58
	% EndTime: 2019-12-29 19:04:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->3), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t16 = qJ(1) + qJ(2);
	t15 = cos(t16);
	t14 = sin(t16);
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t14, t14, 0, 0, 0; -t15, -t15, 0, 0, 0; 0, 0, 0, 0, 0; t15, t15, 0, 0, 0; t14, t14, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:04:51
	% EndTime: 2019-12-29 19:04:51
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->7), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t32 = cos(qJ(4));
	t31 = sin(qJ(4));
	t30 = qJ(1) + qJ(2);
	t29 = cos(t30);
	t28 = sin(t30);
	t27 = t29 * t32;
	t26 = t29 * t31;
	t25 = t28 * t32;
	t24 = t28 * t31;
	t1 = [t26, t26, 0, t25, 0; t24, t24, 0, -t27, 0; 0, 0, 0, -t31, 0; t27, t27, 0, -t24, 0; t25, t25, 0, t26, 0; 0, 0, 0, -t32, 0; -t28, -t28, 0, 0, 0; t29, t29, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:04:58
	% EndTime: 2019-12-29 19:04:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (50->12), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->11)
	t46 = qJ(1) + qJ(2);
	t45 = qJ(4) + qJ(5);
	t44 = cos(t46);
	t43 = cos(t45);
	t42 = sin(t46);
	t41 = sin(t45);
	t40 = t44 * t43;
	t39 = t44 * t41;
	t38 = t42 * t43;
	t37 = t42 * t41;
	t1 = [t39, t39, 0, t38, t38; t37, t37, 0, -t40, -t40; 0, 0, 0, -t41, -t41; t40, t40, 0, -t37, -t37; t38, t38, 0, t39, t39; 0, 0, 0, -t43, -t43; -t42, -t42, 0, 0, 0; t44, t44, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end