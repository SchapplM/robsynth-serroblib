% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR5
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
%   Siehe auch: S5RRPRR5_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRPRR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:03:43
	% EndTime: 2022-01-20 11:03:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:03:43
	% EndTime: 2022-01-20 11:03:43
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
	% StartTime: 2022-01-20 11:03:43
	% EndTime: 2022-01-20 11:03:43
	% DurationCPUTime: 0.02s
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
	% StartTime: 2022-01-20 11:03:43
	% EndTime: 2022-01-20 11:03:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t23 = qJ(1) + qJ(2);
	t21 = sin(t23);
	t25 = cos(pkin(9));
	t27 = t21 * t25;
	t22 = cos(t23);
	t24 = sin(pkin(9));
	t26 = t22 * t24;
	t20 = t22 * t25;
	t19 = t21 * t24;
	t1 = [-t27, -t27, 0, 0, 0; t20, t20, 0, 0, 0; 0, 0, 0, 0, 0; t19, t19, 0, 0, 0; -t26, -t26, 0, 0, 0; 0, 0, 0, 0, 0; t22, t22, 0, 0, 0; t21, t21, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:03:43
	% EndTime: 2022-01-20 11:03:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (39->11), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t29 = pkin(9) + qJ(4);
	t26 = cos(t29);
	t30 = qJ(1) + qJ(2);
	t27 = sin(t30);
	t32 = t27 * t26;
	t25 = sin(t29);
	t28 = cos(t30);
	t31 = t28 * t25;
	t24 = t28 * t26;
	t23 = t27 * t25;
	t1 = [-t32, -t32, 0, -t31, 0; t24, t24, 0, -t23, 0; 0, 0, 0, t26, 0; t23, t23, 0, -t24, 0; -t31, -t31, 0, -t32, 0; 0, 0, 0, -t25, 0; t28, t28, 0, 0, 0; t27, t27, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:03:43
	% EndTime: 2022-01-20 11:03:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (74->16), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->11)
	t35 = pkin(9) + qJ(4) + qJ(5);
	t34 = cos(t35);
	t38 = qJ(1) + qJ(2);
	t36 = sin(t38);
	t40 = t36 * t34;
	t33 = sin(t35);
	t37 = cos(t38);
	t39 = t37 * t33;
	t32 = t37 * t34;
	t31 = t36 * t33;
	t1 = [-t40, -t40, 0, -t39, -t39; t32, t32, 0, -t31, -t31; 0, 0, 0, t34, t34; t31, t31, 0, -t32, -t32; -t39, -t39, 0, -t40, -t40; 0, 0, 0, -t33, -t33; t37, t37, 0, 0, 0; t36, t36, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
end