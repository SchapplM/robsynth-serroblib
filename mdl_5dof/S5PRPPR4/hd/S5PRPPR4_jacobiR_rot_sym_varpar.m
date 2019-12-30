% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPPR4
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRPPR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:28:42
	% EndTime: 2019-12-29 15:28:42
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:28:48
	% EndTime: 2019-12-29 15:28:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:28:43
	% EndTime: 2019-12-29 15:28:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = pkin(7) + qJ(2);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [0, -t10, 0, 0, 0; 0, t11, 0, 0, 0; 0, 0, 0, 0, 0; 0, -t11, 0, 0, 0; 0, -t10, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:28:43
	% EndTime: 2019-12-29 15:28:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->3), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->6)
	t12 = cos(pkin(8));
	t11 = sin(pkin(8));
	t10 = pkin(7) + qJ(2);
	t9 = cos(t10);
	t8 = sin(t10);
	t1 = [0, -t8 * t12, 0, 0, 0; 0, t9 * t12, 0, 0, 0; 0, 0, 0, 0, 0; 0, t8 * t11, 0, 0, 0; 0, -t9 * t11, 0, 0, 0; 0, 0, 0, 0, 0; 0, t9, 0, 0, 0; 0, t8, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:28:43
	% EndTime: 2019-12-29 15:28:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->3), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->6)
	t37 = cos(pkin(8));
	t36 = sin(pkin(8));
	t35 = pkin(7) + qJ(2);
	t34 = cos(t35);
	t33 = sin(t35);
	t1 = [0, -t33 * t37, 0, 0, 0; 0, t34 * t37, 0, 0, 0; 0, 0, 0, 0, 0; 0, t34, 0, 0, 0; 0, t33, 0, 0, 0; 0, 0, 0, 0, 0; 0, -t33 * t36, 0, 0, 0; 0, t34 * t36, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:28:42
	% EndTime: 2019-12-29 15:28:43
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (30->11), mult. (36->8), div. (0->0), fcn. (58->6), ass. (0->14)
	t26 = sin(pkin(8));
	t27 = cos(pkin(8));
	t28 = sin(qJ(5));
	t29 = cos(qJ(5));
	t31 = t26 * t29 - t27 * t28;
	t30 = t26 * t28 + t27 * t29;
	t25 = pkin(7) + qJ(2);
	t24 = cos(t25);
	t23 = sin(t25);
	t22 = t30 * t24;
	t21 = t31 * t24;
	t20 = t30 * t23;
	t19 = t31 * t23;
	t1 = [0, -t20, 0, 0, t21; 0, t22, 0, 0, t19; 0, 0, 0, 0, -t30; 0, -t19, 0, 0, -t22; 0, t21, 0, 0, -t20; 0, 0, 0, 0, -t31; 0, -t24, 0, 0, 0; 0, -t23, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end