% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRP10
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
%   Wie in S5RRPRP10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRPRP10_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP10_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_jacobia_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:11:30
	% EndTime: 2019-12-31 20:11:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:11:30
	% EndTime: 2019-12-31 20:11:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:11:30
	% EndTime: 2019-12-31 20:11:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:11:30
	% EndTime: 2019-12-31 20:11:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (82->16), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->27)
	t29 = cos(qJ(1));
	t25 = t29 ^ 2;
	t28 = cos(qJ(2));
	t26 = sin(qJ(2));
	t27 = sin(qJ(1));
	t32 = t27 * t26;
	t18 = atan2(-t32, -t28);
	t16 = sin(t18);
	t17 = cos(t18);
	t14 = -t16 * t32 - t17 * t28;
	t13 = 0.1e1 / t14 ^ 2;
	t38 = t13 * t26;
	t37 = t16 * t28;
	t21 = t26 ^ 2;
	t31 = t28 ^ 2;
	t36 = t21 / t31;
	t30 = t27 ^ 2;
	t35 = 0.1e1 / t30 * t25;
	t34 = t26 * t29;
	t19 = 0.1e1 / (t30 * t36 + 0.1e1);
	t33 = t27 * t19;
	t23 = 0.1e1 / t28;
	t20 = 0.1e1 / (t31 * t35 + 0.1e1);
	t15 = (0.1e1 + t36) * t33;
	t12 = 0.1e1 / t14;
	t11 = 0.1e1 / (t25 * t21 * t13 + 0.1e1);
	t1 = [t23 * t19 * t34, t15, 0, 0, 0; (-t12 * t32 - (-t17 * t21 * t23 * t33 + (t19 - 0.1e1) * t26 * t16) * t25 * t38) * t11, (t28 * t12 - (-t27 * t37 + t17 * t26 + (-t17 * t32 + t37) * t15) * t38) * t29 * t11, 0, 0, 0; (-0.1e1 - t35) * t28 * t20, -0.1e1 / t27 * t20 * t34, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:11:30
	% EndTime: 2019-12-31 20:11:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (87->20), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->35)
	t40 = sin(qJ(2));
	t41 = sin(qJ(1));
	t43 = cos(qJ(2));
	t49 = t41 * t43;
	t34 = atan2(-t49, t40);
	t32 = sin(t34);
	t33 = cos(t34);
	t25 = -t32 * t49 + t33 * t40;
	t24 = 0.1e1 / t25 ^ 2;
	t44 = cos(qJ(1));
	t56 = t24 * t44 ^ 2;
	t39 = sin(qJ(4));
	t47 = t44 * t39;
	t42 = cos(qJ(4));
	t50 = t41 * t42;
	t31 = t40 * t47 + t50;
	t29 = 0.1e1 / t31 ^ 2;
	t46 = t44 * t42;
	t51 = t41 * t39;
	t30 = -t40 * t46 + t51;
	t55 = t29 * t30;
	t54 = t32 * t40;
	t38 = t43 ^ 2;
	t53 = 0.1e1 / t40 ^ 2 * t38;
	t35 = 0.1e1 / (t41 ^ 2 * t53 + 0.1e1);
	t52 = t41 * t35;
	t48 = t43 * t44;
	t45 = t30 ^ 2 * t29 + 0.1e1;
	t36 = 0.1e1 / t40;
	t28 = 0.1e1 / t31;
	t27 = (0.1e1 + t53) * t52;
	t26 = 0.1e1 / t45;
	t23 = 0.1e1 / t25;
	t22 = 0.1e1 / (t38 * t56 + 0.1e1);
	t1 = [-t36 * t35 * t48, t27, 0, 0, 0; (-t23 * t49 - (t33 * t36 * t38 * t52 + (t35 - 0.1e1) * t43 * t32) * t43 * t56) * t22, (-t40 * t23 - (t41 * t54 + t33 * t43 + (-t33 * t49 - t54) * t27) * t43 * t24) * t44 * t22, 0, 0, 0; ((t40 * t50 + t47) * t28 - (-t40 * t51 + t46) * t55) * t26, (-t28 * t42 - t39 * t55) * t26 * t48, 0, t45 * t26, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:11:31
	% EndTime: 2019-12-31 20:11:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (87->20), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->35)
	t44 = sin(qJ(2));
	t45 = sin(qJ(1));
	t47 = cos(qJ(2));
	t53 = t45 * t47;
	t38 = atan2(-t53, t44);
	t36 = sin(t38);
	t37 = cos(t38);
	t29 = -t36 * t53 + t37 * t44;
	t28 = 0.1e1 / t29 ^ 2;
	t48 = cos(qJ(1));
	t60 = t28 * t48 ^ 2;
	t43 = sin(qJ(4));
	t51 = t48 * t43;
	t46 = cos(qJ(4));
	t54 = t45 * t46;
	t35 = t44 * t51 + t54;
	t33 = 0.1e1 / t35 ^ 2;
	t50 = t48 * t46;
	t55 = t45 * t43;
	t34 = -t44 * t50 + t55;
	t59 = t33 * t34;
	t58 = t36 * t44;
	t42 = t47 ^ 2;
	t57 = 0.1e1 / t44 ^ 2 * t42;
	t39 = 0.1e1 / (t45 ^ 2 * t57 + 0.1e1);
	t56 = t45 * t39;
	t52 = t47 * t48;
	t49 = t34 ^ 2 * t33 + 0.1e1;
	t40 = 0.1e1 / t44;
	t32 = 0.1e1 / t35;
	t31 = (0.1e1 + t57) * t56;
	t30 = 0.1e1 / t49;
	t27 = 0.1e1 / t29;
	t26 = 0.1e1 / (t42 * t60 + 0.1e1);
	t1 = [-t40 * t39 * t52, t31, 0, 0, 0; (-t27 * t53 - (t37 * t40 * t42 * t56 + (t39 - 0.1e1) * t47 * t36) * t47 * t60) * t26, (-t44 * t27 - (t45 * t58 + t37 * t47 + (-t37 * t53 - t58) * t31) * t47 * t28) * t48 * t26, 0, 0, 0; ((t44 * t54 + t51) * t32 - (-t44 * t55 + t50) * t59) * t30, (-t32 * t46 - t43 * t59) * t30 * t52, 0, t49 * t30, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end