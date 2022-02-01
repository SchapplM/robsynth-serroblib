% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR5
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
%   Siehe auch: S5RRRRR5_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:45
	% EndTime: 2022-01-20 12:02:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:45
	% EndTime: 2022-01-20 12:02:45
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
	% StartTime: 2022-01-20 12:02:45
	% EndTime: 2022-01-20 12:02:45
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
	% StartTime: 2022-01-20 12:02:45
	% EndTime: 2022-01-20 12:02:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->10), mult. (0->0), div. (0->0), fcn. (12->2), ass. (0->4)
	t21 = qJ(1) + qJ(2) + qJ(3);
	t20 = cos(t21);
	t19 = sin(t21);
	t1 = [-t19, -t19, -t19, 0, 0; t20, t20, t20, 0, 0; 0, 0, 0, 0, 0; -t20, -t20, -t20, 0, 0; -t19, -t19, -t19, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:45
	% EndTime: 2022-01-20 12:02:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (55->12), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t31 = qJ(1) + qJ(2) + qJ(3);
	t29 = sin(t31);
	t33 = cos(qJ(4));
	t35 = t29 * t33;
	t30 = cos(t31);
	t32 = sin(qJ(4));
	t34 = t30 * t32;
	t28 = t30 * t33;
	t27 = t29 * t32;
	t1 = [-t35, -t35, -t35, -t34, 0; t28, t28, t28, -t27, 0; 0, 0, 0, t33, 0; t27, t27, t27, -t28, 0; -t34, -t34, -t34, -t35, 0; 0, 0, 0, -t32, 0; t30, t30, t30, 0, 0; t29, t29, t29, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 12:02:45
	% EndTime: 2022-01-20 12:02:45
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (92->18), mult. (20->4), div. (0->0), fcn. (50->4), ass. (0->11)
	t40 = qJ(1) + qJ(2) + qJ(3);
	t36 = sin(t40);
	t41 = qJ(4) + qJ(5);
	t39 = cos(t41);
	t43 = t36 * t39;
	t37 = cos(t40);
	t38 = sin(t41);
	t42 = t37 * t38;
	t35 = t37 * t39;
	t34 = t36 * t38;
	t1 = [-t43, -t43, -t43, -t42, -t42; t35, t35, t35, -t34, -t34; 0, 0, 0, t39, t39; t34, t34, t34, -t35, -t35; -t42, -t42, -t42, -t43, -t43; 0, 0, 0, -t38, -t38; t37, t37, t37, 0, 0; t36, t36, t36, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
end