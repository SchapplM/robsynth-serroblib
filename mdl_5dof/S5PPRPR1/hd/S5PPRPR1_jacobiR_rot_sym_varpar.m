% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRPR1
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PPRPR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->6), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->6)
	t11 = cos(pkin(7));
	t10 = sin(pkin(7));
	t9 = pkin(8) + qJ(3);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, -t11 * t7, 0, 0; 0, 0, -t10 * t7, 0, 0; 0, 0, t8, 0, 0; 0, 0, -t11 * t8, 0, 0; 0, 0, -t10 * t8, 0, 0; 0, 0, -t7, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->4), mult. (12->10), div. (0->0), fcn. (21->6), ass. (0->10)
	t38 = pkin(8) + qJ(3);
	t36 = sin(t38);
	t40 = sin(pkin(7));
	t44 = t36 * t40;
	t42 = cos(pkin(7));
	t43 = t36 * t42;
	t41 = cos(pkin(9));
	t39 = sin(pkin(9));
	t37 = cos(t38);
	t1 = [0, 0, -t41 * t43, 0, 0; 0, 0, -t41 * t44, 0, 0; 0, 0, t37 * t41, 0, 0; 0, 0, t39 * t43, 0, 0; 0, 0, t39 * t44, 0, 0; 0, 0, -t37 * t39, 0, 0; 0, 0, t42 * t37, 0, 0; 0, 0, t40 * t37, 0, 0; 0, 0, t36, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (40->11), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->15)
	t53 = pkin(9) + qJ(5);
	t49 = sin(t53);
	t55 = sin(pkin(7));
	t62 = t55 * t49;
	t51 = cos(t53);
	t61 = t55 * t51;
	t54 = pkin(8) + qJ(3);
	t52 = cos(t54);
	t60 = t55 * t52;
	t56 = cos(pkin(7));
	t59 = t56 * t49;
	t58 = t56 * t51;
	t57 = t56 * t52;
	t50 = sin(t54);
	t1 = [0, 0, -t50 * t58, 0, -t49 * t57 + t61; 0, 0, -t50 * t61, 0, -t49 * t60 - t58; 0, 0, t52 * t51, 0, -t50 * t49; 0, 0, t50 * t59, 0, -t51 * t57 - t62; 0, 0, t50 * t62, 0, -t51 * t60 + t59; 0, 0, -t52 * t49, 0, -t50 * t51; 0, 0, t57, 0, 0; 0, 0, t60, 0, 0; 0, 0, t50, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end