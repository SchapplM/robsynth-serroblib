% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RPPRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPPRP5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_jacobia_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:09
	% EndTime: 2019-12-31 17:54:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:09
	% EndTime: 2019-12-31 17:54:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:09
	% EndTime: 2019-12-31 17:54:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:09
	% EndTime: 2019-12-31 17:54:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (37->10), mult. (63->23), div. (23->8), fcn. (105->7), ass. (0->19)
	t23 = cos(pkin(7));
	t22 = sin(pkin(7));
	t24 = sin(qJ(1));
	t27 = t24 * t22;
	t14 = atan2(-t27, -t23);
	t12 = sin(t14);
	t13 = cos(t14);
	t11 = -t12 * t27 - t13 * t23;
	t25 = cos(qJ(1));
	t26 = t25 ^ 2;
	t30 = 0.1e1 / t11 ^ 2 * t26;
	t17 = t22 ^ 2;
	t19 = 0.1e1 / t23 ^ 2;
	t20 = t24 ^ 2;
	t15 = 0.1e1 / (t20 * t17 * t19 + 0.1e1);
	t18 = 0.1e1 / t23;
	t29 = t15 * t18;
	t28 = t20 / t26;
	t1 = [t25 * t22 * t29, 0, 0, 0, 0; (-0.1e1 / t11 * t27 - (-t13 * t17 * t24 * t29 + (t15 - 0.1e1) * t22 * t12) * t22 * t30) / (t17 * t30 + 0.1e1), 0, 0, 0, 0; (-0.1e1 - t28) * t18 / (t19 * t28 + 0.1e1), 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:09
	% EndTime: 2019-12-31 17:54:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:09
	% EndTime: 2019-12-31 17:54:10
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (158->16), mult. (454->43), div. (47->9), fcn. (678->9), ass. (0->29)
	t47 = cos(pkin(7));
	t48 = sin(qJ(4));
	t52 = sin(pkin(7));
	t57 = cos(qJ(4));
	t44 = -t47 * t48 + t52 * t57;
	t49 = sin(qJ(1));
	t36 = t44 * t49;
	t43 = t47 * t57 + t52 * t48;
	t33 = atan2(t36, t43);
	t31 = cos(t33);
	t56 = t31 * t36;
	t50 = cos(qJ(1));
	t40 = t43 * t50;
	t35 = 0.1e1 / t40 ^ 2;
	t55 = t35 * t49;
	t30 = sin(t33);
	t29 = t30 * t36 + t31 * t43;
	t28 = 0.1e1 / t29 ^ 2;
	t39 = t44 * t50;
	t54 = t39 ^ 2 * t28;
	t42 = 0.1e1 / t43 ^ 2;
	t41 = 0.1e1 / t43;
	t38 = t43 * t49;
	t34 = 0.1e1 / (t49 ^ 2 * t35 + 0.1e1);
	t32 = 0.1e1 / (t36 ^ 2 * t42 + 0.1e1);
	t27 = 0.1e1 / t29;
	t26 = 0.1e1 / (0.1e1 + t54);
	t25 = (-t36 * t42 * t44 - t38 * t41) * t32;
	t1 = [t39 * t41 * t32, 0, 0, t25, 0; (t36 * t27 - (-t30 + (-t41 * t56 + t30) * t32) * t54) * t26, 0, 0, (t40 * t27 + (-t30 * t38 + t31 * t44 + (-t30 * t43 + t56) * t25) * t39 * t28) * t26, 0; (t50 / t40 + t38 * t55) * t34, 0, 0, -t39 * t34 * t55, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end