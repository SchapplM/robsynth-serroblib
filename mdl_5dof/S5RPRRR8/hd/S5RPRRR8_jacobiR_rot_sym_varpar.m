% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 17:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRRR8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:48:55
	% EndTime: 2019-12-29 17:48:55
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:48:55
	% EndTime: 2019-12-29 17:48:55
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
	% StartTime: 2019-12-29 17:48:56
	% EndTime: 2019-12-29 17:48:56
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
	% StartTime: 2019-12-29 17:48:55
	% EndTime: 2019-12-29 17:48:55
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (8->6), mult. (16->4), div. (0->0), fcn. (32->4), ass. (0->7)
	t35 = cos(qJ(1));
	t34 = cos(qJ(3));
	t33 = sin(qJ(1));
	t32 = sin(qJ(3));
	t28 = t35 * t32 - t33 * t34;
	t27 = -t33 * t32 - t35 * t34;
	t1 = [t28, 0, -t28, 0, 0; -t27, 0, t27, 0, 0; 0, 0, 0, 0, 0; -t27, 0, t27, 0, 0; -t28, 0, t28, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:49:00
	% EndTime: 2019-12-29 17:49:01
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (25->9), mult. (44->8), div. (0->0), fcn. (78->6), ass. (0->13)
	t37 = cos(qJ(1));
	t36 = cos(qJ(3));
	t35 = sin(qJ(1));
	t34 = sin(qJ(3));
	t23 = -t35 * t34 - t37 * t36;
	t28 = sin(qJ(4));
	t33 = t23 * t28;
	t29 = cos(qJ(4));
	t32 = t23 * t29;
	t24 = t37 * t34 - t35 * t36;
	t31 = t24 * t28;
	t30 = t24 * t29;
	t1 = [t30, 0, -t30, t33, 0; -t32, 0, t32, t31, 0; 0, 0, 0, -t29, 0; -t31, 0, t31, t32, 0; t33, 0, -t33, t30, 0; 0, 0, 0, t28, 0; t23, 0, -t23, 0, 0; t24, 0, -t24, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:48:55
	% EndTime: 2019-12-29 17:48:55
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (54->11), mult. (56->8), div. (0->0), fcn. (100->6), ass. (0->14)
	t50 = cos(qJ(1));
	t49 = cos(qJ(3));
	t48 = sin(qJ(1));
	t47 = sin(qJ(3));
	t35 = -t48 * t47 - t50 * t49;
	t42 = qJ(4) + qJ(5);
	t40 = sin(t42);
	t46 = t35 * t40;
	t41 = cos(t42);
	t45 = t35 * t41;
	t36 = t50 * t47 - t48 * t49;
	t44 = t36 * t40;
	t43 = t36 * t41;
	t1 = [t43, 0, -t43, t46, t46; -t45, 0, t45, t44, t44; 0, 0, 0, -t41, -t41; -t44, 0, t44, t45, t45; t46, 0, -t46, t43, t43; 0, 0, 0, t40, t40; t35, 0, -t35, 0, 0; t36, 0, -t36, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end