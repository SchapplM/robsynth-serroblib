% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRPP3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPP3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_jacobiR_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:38:03
	% EndTime: 2019-12-29 16:38:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:38:03
	% EndTime: 2019-12-29 16:38:03
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
	% StartTime: 2019-12-29 16:38:03
	% EndTime: 2019-12-29 16:38:03
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t10 = cos(qJ(1));
	t9 = sin(qJ(1));
	t8 = cos(pkin(7));
	t7 = sin(pkin(7));
	t1 = [-t9 * t8, 0, 0, 0, 0; t10 * t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9 * t7, 0, 0, 0, 0; -t10 * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0; t10, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:38:03
	% EndTime: 2019-12-29 16:38:03
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t14 = pkin(7) + qJ(3);
	t12 = sin(t14);
	t15 = sin(qJ(1));
	t20 = t15 * t12;
	t13 = cos(t14);
	t19 = t15 * t13;
	t16 = cos(qJ(1));
	t18 = t16 * t12;
	t17 = t16 * t13;
	t1 = [-t19, 0, -t18, 0, 0; t17, 0, -t20, 0, 0; 0, 0, t13, 0, 0; t20, 0, -t17, 0, 0; -t18, 0, -t19, 0, 0; 0, 0, -t12, 0, 0; t16, 0, 0, 0, 0; t15, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:38:09
	% EndTime: 2019-12-29 16:38:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (13->4), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t47 = pkin(7) + qJ(3);
	t45 = sin(t47);
	t48 = sin(qJ(1));
	t51 = t48 * t45;
	t46 = cos(t47);
	t49 = cos(qJ(1));
	t50 = t49 * t46;
	t44 = t49 * t45;
	t43 = t48 * t46;
	t1 = [t49, 0, 0, 0, 0; t48, 0, 0, 0, 0; 0, 0, 0, 0, 0; t43, 0, t44, 0, 0; -t50, 0, t51, 0, 0; 0, 0, -t46, 0, 0; -t51, 0, t50, 0, 0; t44, 0, t43, 0, 0; 0, 0, t45, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:38:03
	% EndTime: 2019-12-29 16:38:03
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (14->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t46 = pkin(7) + qJ(3);
	t44 = sin(t46);
	t47 = sin(qJ(1));
	t51 = t47 * t44;
	t45 = cos(t46);
	t50 = t47 * t45;
	t48 = cos(qJ(1));
	t49 = t48 * t44;
	t43 = t48 * t45;
	t1 = [t48, 0, 0, 0, 0; t47, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t51, 0, t43, 0, 0; t49, 0, t50, 0, 0; 0, 0, t44, 0, 0; -t50, 0, -t49, 0, 0; t43, 0, -t51, 0, 0; 0, 0, t45, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end