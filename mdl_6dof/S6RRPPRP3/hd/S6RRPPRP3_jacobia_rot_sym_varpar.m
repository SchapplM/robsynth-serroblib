% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RRPPRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:30
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:30:29
	% EndTime: 2019-10-10 09:30:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:30:29
	% EndTime: 2019-10-10 09:30:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:30:29
	% EndTime: 2019-10-10 09:30:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:30:29
	% EndTime: 2019-10-10 09:30:29
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (82->16), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->26)
	t29 = cos(qJ(1));
	t30 = t29 ^ 2;
	t28 = cos(qJ(2));
	t26 = sin(qJ(2));
	t27 = sin(qJ(1));
	t31 = t27 * t26;
	t18 = atan2(-t31, -t28);
	t16 = sin(t18);
	t17 = cos(t18);
	t14 = -t16 * t31 - t17 * t28;
	t13 = 0.1e1 / t14 ^ 2;
	t36 = t13 * t26;
	t35 = t16 * t28;
	t21 = t26 ^ 2;
	t24 = 0.1e1 / t28 ^ 2;
	t34 = t21 * t24;
	t22 = t27 ^ 2;
	t33 = t22 / t30;
	t19 = 0.1e1 / (t22 * t34 + 0.1e1);
	t32 = t27 * t19;
	t23 = 0.1e1 / t28;
	t20 = 0.1e1 / (t24 * t33 + 0.1e1);
	t15 = (0.1e1 + t34) * t32;
	t12 = 0.1e1 / t14;
	t11 = 0.1e1 / (t30 * t21 * t13 + 0.1e1);
	t1 = [t29 * t26 * t23 * t19, t15, 0, 0, 0, 0; (-t12 * t31 - (-t17 * t21 * t23 * t32 + (t19 - 0.1e1) * t26 * t16) * t30 * t36) * t11, (t28 * t12 - (-t27 * t35 + t17 * t26 + (-t17 * t31 + t35) * t15) * t36) * t29 * t11, 0, 0, 0, 0; (-0.1e1 - t33) * t23 * t20, -0.1e1 / t29 * t24 * t20 * t31, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:30:29
	% EndTime: 2019-10-10 09:30:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:30:29
	% EndTime: 2019-10-10 09:30:29
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (87->20), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->35)
	t41 = sin(qJ(2));
	t42 = sin(qJ(1));
	t44 = cos(qJ(2));
	t50 = t42 * t44;
	t35 = atan2(-t50, t41);
	t33 = sin(t35);
	t34 = cos(t35);
	t26 = -t33 * t50 + t34 * t41;
	t25 = 0.1e1 / t26 ^ 2;
	t45 = cos(qJ(1));
	t57 = t25 * t45 ^ 2;
	t43 = cos(qJ(5));
	t47 = t45 * t43;
	t40 = sin(qJ(5));
	t52 = t42 * t40;
	t32 = t41 * t47 - t52;
	t30 = 0.1e1 / t32 ^ 2;
	t48 = t45 * t40;
	t51 = t42 * t43;
	t31 = t41 * t48 + t51;
	t56 = t30 * t31;
	t55 = t33 * t41;
	t39 = t44 ^ 2;
	t54 = 0.1e1 / t41 ^ 2 * t39;
	t36 = 0.1e1 / (t42 ^ 2 * t54 + 0.1e1);
	t53 = t42 * t36;
	t49 = t44 * t45;
	t46 = t31 ^ 2 * t30 + 0.1e1;
	t37 = 0.1e1 / t41;
	t29 = 0.1e1 / t32;
	t28 = (0.1e1 + t54) * t53;
	t27 = 0.1e1 / t46;
	t24 = 0.1e1 / t26;
	t23 = 0.1e1 / (t39 * t57 + 0.1e1);
	t1 = [-t37 * t36 * t49, t28, 0, 0, 0, 0; (-t24 * t50 - (t34 * t37 * t39 * t53 + (t36 - 0.1e1) * t44 * t33) * t44 * t57) * t23, (-t41 * t24 - (t42 * t55 + t34 * t44 + (-t34 * t50 - t55) * t28) * t44 * t25) * t45 * t23, 0, 0, 0, 0; ((-t41 * t52 + t47) * t29 - (-t41 * t51 - t48) * t56) * t27, (t29 * t40 - t43 * t56) * t27 * t49, 0, 0, t46 * t27, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:30:29
	% EndTime: 2019-10-10 09:30:29
	% DurationCPUTime: 0.10s
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
	t46 = cos(qJ(5));
	t50 = t48 * t46;
	t43 = sin(qJ(5));
	t55 = t45 * t43;
	t35 = t44 * t50 - t55;
	t33 = 0.1e1 / t35 ^ 2;
	t51 = t48 * t43;
	t54 = t45 * t46;
	t34 = t44 * t51 + t54;
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
	t1 = [-t40 * t39 * t52, t31, 0, 0, 0, 0; (-t27 * t53 - (t37 * t40 * t42 * t56 + (t39 - 0.1e1) * t47 * t36) * t47 * t60) * t26, (-t44 * t27 - (t45 * t58 + t37 * t47 + (-t37 * t53 - t58) * t31) * t47 * t28) * t48 * t26, 0, 0, 0, 0; ((-t44 * t55 + t50) * t32 - (-t44 * t54 - t51) * t59) * t30, (t32 * t43 - t46 * t59) * t30 * t52, 0, 0, t49 * t30, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end