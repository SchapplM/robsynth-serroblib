% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRPP5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPP5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_jacobiR_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:43:30
	% EndTime: 2019-12-29 16:43:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:43:30
	% EndTime: 2019-12-29 16:43:30
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
	% StartTime: 2019-12-29 16:43:30
	% EndTime: 2019-12-29 16:43:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t5, 0, 0, 0, 0; -t6, 0, 0, 0, 0; 0, 0, 0, 0, 0; t6, 0, 0, 0, 0; t5, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:43:30
	% EndTime: 2019-12-29 16:43:30
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t11 = sin(qJ(3));
	t12 = sin(qJ(1));
	t16 = t12 * t11;
	t13 = cos(qJ(3));
	t14 = cos(qJ(1));
	t15 = t14 * t13;
	t10 = t14 * t11;
	t9 = t12 * t13;
	t1 = [t10, 0, t9, 0, 0; t16, 0, -t15, 0, 0; 0, 0, -t11, 0, 0; t15, 0, -t16, 0, 0; t9, 0, t10, 0, 0; 0, 0, -t13, 0, 0; -t12, 0, 0, 0, 0; t14, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:43:35
	% EndTime: 2019-12-29 16:43:35
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t39 = sin(qJ(1));
	t40 = cos(qJ(3));
	t44 = t39 * t40;
	t38 = sin(qJ(3));
	t41 = cos(qJ(1));
	t43 = t41 * t38;
	t42 = t41 * t40;
	t37 = t39 * t38;
	t1 = [t43, 0, t44, 0, 0; t37, 0, -t42, 0, 0; 0, 0, -t38, 0, 0; -t39, 0, 0, 0, 0; t41, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t42, 0, t37, 0, 0; -t44, 0, -t43, 0, 0; 0, 0, t40, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:43:35
	% EndTime: 2019-12-29 16:43:35
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t13 = sin(qJ(1));
	t14 = cos(qJ(3));
	t18 = t13 * t14;
	t12 = sin(qJ(3));
	t15 = cos(qJ(1));
	t17 = t15 * t12;
	t16 = t15 * t14;
	t11 = t13 * t12;
	t1 = [t17, 0, t18, 0, 0; t11, 0, -t16, 0, 0; 0, 0, -t12, 0, 0; -t16, 0, t11, 0, 0; -t18, 0, -t17, 0, 0; 0, 0, t14, 0, 0; t13, 0, 0, 0, 0; -t15, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end