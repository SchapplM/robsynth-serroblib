% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRRR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:15:07
	% EndTime: 2019-12-05 18:15:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:15:06
	% EndTime: 2019-12-05 18:15:06
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
	% StartTime: 2019-12-05 18:15:06
	% EndTime: 2019-12-05 18:15:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t6 = qJ(1) + pkin(9);
	t5 = cos(t6);
	t4 = sin(t6);
	t1 = [0, 0, 0, 0, 0; -t5, 0, 0, 0, 0; -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; t4, 0, 0, 0, 0; -t5, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:15:06
	% EndTime: 2019-12-05 18:15:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t10 = qJ(1) + pkin(9) + qJ(3);
	t9 = cos(t10);
	t8 = sin(t10);
	t1 = [0, 0, 0, 0, 0; -t9, 0, -t9, 0, 0; -t8, 0, -t8, 0, 0; 0, 0, 0, 0, 0; t8, 0, t8, 0, 0; -t9, 0, -t9, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:15:06
	% EndTime: 2019-12-05 18:15:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (45->10), mult. (0->0), div. (0->0), fcn. (12->2), ass. (0->4)
	t15 = qJ(1) + pkin(9) + qJ(3) + qJ(4);
	t14 = cos(t15);
	t13 = sin(t15);
	t1 = [0, 0, 0, 0, 0; -t14, 0, -t14, -t14, 0; -t13, 0, -t13, -t13, 0; 0, 0, 0, 0, 0; t13, 0, t13, t13, 0; -t14, 0, -t14, -t14, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:15:06
	% EndTime: 2019-12-05 18:15:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (78->13), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t42 = qJ(1) + pkin(9) + qJ(3) + qJ(4);
	t40 = sin(t42);
	t44 = cos(qJ(5));
	t46 = t40 * t44;
	t41 = cos(t42);
	t45 = t41 * t44;
	t43 = sin(qJ(5));
	t39 = t41 * t43;
	t38 = t40 * t43;
	t1 = [0, 0, 0, 0, t44; -t45, 0, -t45, -t45, t38; -t46, 0, -t46, -t46, -t39; 0, 0, 0, 0, -t43; t39, 0, t39, t39, t46; t38, 0, t38, t38, -t45; 0, 0, 0, 0, 0; -t40, 0, -t40, -t40, 0; t41, 0, t41, t41, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end