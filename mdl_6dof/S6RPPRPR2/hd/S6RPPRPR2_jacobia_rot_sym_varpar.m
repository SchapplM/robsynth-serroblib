% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR2
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
%   Wie in S6RPPRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:36
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (263->18), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->29)
	t35 = qJ(1) + pkin(9);
	t33 = cos(t35);
	t29 = t33 ^ 2;
	t34 = pkin(10) + qJ(4);
	t32 = cos(t34);
	t30 = sin(t34);
	t31 = sin(t35);
	t38 = t31 * t30;
	t22 = atan2(-t38, -t32);
	t20 = sin(t22);
	t21 = cos(t22);
	t18 = -t20 * t38 - t21 * t32;
	t17 = 0.1e1 / t18 ^ 2;
	t44 = t17 * t30;
	t43 = t20 * t32;
	t25 = t30 ^ 2;
	t37 = t32 ^ 2;
	t42 = t25 / t37;
	t36 = t31 ^ 2;
	t41 = 0.1e1 / t36 * t29;
	t40 = t30 * t33;
	t23 = 0.1e1 / (t36 * t42 + 0.1e1);
	t39 = t31 * t23;
	t27 = 0.1e1 / t32;
	t24 = 0.1e1 / (t37 * t41 + 0.1e1);
	t19 = (0.1e1 + t42) * t39;
	t16 = 0.1e1 / t18;
	t15 = 0.1e1 / (t29 * t25 * t17 + 0.1e1);
	t1 = [t27 * t23 * t40, 0, 0, t19, 0, 0; (-t16 * t38 - (-t21 * t25 * t27 * t39 + (t23 - 0.1e1) * t30 * t20) * t29 * t44) * t15, 0, 0, (t32 * t16 - (-t31 * t43 + t21 * t30 + (-t21 * t38 + t43) * t19) * t44) * t33 * t15, 0, 0; (-0.1e1 - t41) * t32 * t24, 0, 0, -0.1e1 / t31 * t24 * t40, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:36:00
	% EndTime: 2019-10-09 23:36:00
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (325->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t47 = pkin(10) + qJ(4);
	t43 = sin(t47);
	t48 = qJ(1) + pkin(9);
	t44 = sin(t48);
	t45 = cos(t47);
	t57 = t44 * t45;
	t38 = atan2(-t57, t43);
	t36 = sin(t38);
	t37 = cos(t38);
	t30 = -t36 * t57 + t37 * t43;
	t28 = 0.1e1 / t30 ^ 2;
	t46 = cos(t48);
	t62 = t28 * t46 ^ 2;
	t49 = sin(qJ(6));
	t53 = t46 * t49;
	t50 = cos(qJ(6));
	t55 = t44 * t50;
	t35 = t43 * t53 + t55;
	t33 = 0.1e1 / t35 ^ 2;
	t52 = t46 * t50;
	t56 = t44 * t49;
	t34 = -t43 * t52 + t56;
	t61 = t33 * t34;
	t60 = t36 * t43;
	t42 = t45 ^ 2;
	t59 = 0.1e1 / t43 ^ 2 * t42;
	t39 = 0.1e1 / (t44 ^ 2 * t59 + 0.1e1);
	t58 = t44 * t39;
	t54 = t45 * t46;
	t51 = t34 ^ 2 * t33 + 0.1e1;
	t40 = 0.1e1 / t43;
	t32 = 0.1e1 / t35;
	t31 = 0.1e1 / t51;
	t29 = (0.1e1 + t59) * t58;
	t27 = 0.1e1 / t30;
	t26 = 0.1e1 / (t42 * t62 + 0.1e1);
	t1 = [-t40 * t39 * t54, 0, 0, t29, 0, 0; (-t27 * t57 - (t37 * t40 * t42 * t58 + (t39 - 0.1e1) * t45 * t36) * t45 * t62) * t26, 0, 0, (-t43 * t27 - (t44 * t60 + t37 * t45 + (-t37 * t57 - t60) * t29) * t45 * t28) * t46 * t26, 0, 0; ((t43 * t55 + t53) * t32 - (-t43 * t56 + t52) * t61) * t31, 0, 0, (-t32 * t50 - t49 * t61) * t31 * t54, 0, t51 * t31;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end