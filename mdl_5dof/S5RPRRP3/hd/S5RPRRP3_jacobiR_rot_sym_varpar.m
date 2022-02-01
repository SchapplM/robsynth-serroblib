% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRP3
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
%   Siehe auch: S5RPRRP3_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRRP3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(8);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [-t10, 0, 0, 0, 0; t11, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t11, 0, 0, 0, 0; -t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t13 = qJ(1) + pkin(8);
	t11 = sin(t13);
	t14 = sin(qJ(3));
	t19 = t11 * t14;
	t15 = cos(qJ(3));
	t18 = t11 * t15;
	t12 = cos(t13);
	t17 = t12 * t14;
	t16 = t12 * t15;
	t1 = [-t18, 0, -t17, 0, 0; t16, 0, -t19, 0, 0; 0, 0, t15, 0, 0; t19, 0, -t16, 0, 0; -t17, 0, -t18, 0, 0; 0, 0, -t14, 0, 0; t12, 0, 0, 0, 0; t11, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (42->14), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t26 = qJ(1) + pkin(8);
	t22 = sin(t26);
	t27 = qJ(3) + qJ(4);
	t24 = sin(t27);
	t31 = t22 * t24;
	t25 = cos(t27);
	t30 = t22 * t25;
	t23 = cos(t26);
	t29 = t23 * t24;
	t28 = t23 * t25;
	t1 = [-t30, 0, -t29, -t29, 0; t28, 0, -t31, -t31, 0; 0, 0, t25, t25, 0; t31, 0, -t28, -t28, 0; -t29, 0, -t30, -t30, 0; 0, 0, -t24, -t24, 0; t23, 0, 0, 0, 0; t22, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (42->14), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t29 = qJ(1) + pkin(8);
	t25 = sin(t29);
	t30 = qJ(3) + qJ(4);
	t27 = sin(t30);
	t34 = t25 * t27;
	t28 = cos(t30);
	t33 = t25 * t28;
	t26 = cos(t29);
	t32 = t26 * t27;
	t31 = t26 * t28;
	t1 = [-t33, 0, -t32, -t32, 0; t31, 0, -t34, -t34, 0; 0, 0, t28, t28, 0; t34, 0, -t31, -t31, 0; -t32, 0, -t33, -t33, 0; 0, 0, -t27, -t27, 0; t26, 0, 0, 0, 0; t25, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
end