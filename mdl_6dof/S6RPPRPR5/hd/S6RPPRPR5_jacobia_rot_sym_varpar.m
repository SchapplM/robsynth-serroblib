% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR5
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
%   Wie in S6RPPRPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (53->18), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->34)
	t33 = sin(qJ(4));
	t34 = sin(qJ(1));
	t35 = cos(qJ(4));
	t40 = t34 * t35;
	t27 = atan2(t40, t33);
	t24 = sin(t27);
	t25 = cos(t27);
	t18 = t24 * t40 + t25 * t33;
	t17 = 0.1e1 / t18 ^ 2;
	t36 = cos(qJ(1));
	t47 = t17 * t36 ^ 2;
	t32 = cos(pkin(9));
	t37 = t36 * t32;
	t31 = sin(pkin(9));
	t42 = t34 * t31;
	t23 = t33 * t37 - t42;
	t21 = 0.1e1 / t23 ^ 2;
	t38 = t36 * t31;
	t41 = t34 * t32;
	t22 = t33 * t38 + t41;
	t46 = t21 * t22;
	t45 = t24 * t33;
	t30 = t35 ^ 2;
	t44 = 0.1e1 / t33 ^ 2 * t30;
	t26 = 0.1e1 / (t34 ^ 2 * t44 + 0.1e1);
	t43 = t34 * t26;
	t39 = t35 * t36;
	t28 = 0.1e1 / t33;
	t20 = 0.1e1 / t23;
	t19 = (-0.1e1 - t44) * t43;
	t16 = 0.1e1 / t18;
	t15 = 0.1e1 / (t22 ^ 2 * t21 + 0.1e1);
	t14 = 0.1e1 / (t30 * t47 + 0.1e1);
	t1 = [t28 * t26 * t39, 0, 0, t19, 0, 0; (t16 * t40 + (t25 * t28 * t30 * t43 + (-t26 + 0.1e1) * t35 * t24) * t35 * t47) * t14, 0, 0, (t33 * t16 + (-t34 * t45 + t25 * t35 + (t25 * t40 - t45) * t19) * t35 * t17) * t36 * t14, 0, 0; ((-t33 * t42 + t37) * t20 - (-t33 * t41 - t38) * t46) * t15, 0, 0, (t20 * t31 - t32 * t46) * t15 * t39, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:41:04
	% EndTime: 2019-10-09 23:41:04
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (111->19), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
	t44 = sin(qJ(4));
	t45 = sin(qJ(1));
	t46 = cos(qJ(4));
	t52 = t45 * t46;
	t37 = atan2(t52, t44);
	t34 = sin(t37);
	t35 = cos(t37);
	t28 = t34 * t52 + t35 * t44;
	t27 = 0.1e1 / t28 ^ 2;
	t47 = cos(qJ(1));
	t59 = t27 * t47 ^ 2;
	t40 = pkin(9) + qJ(6);
	t39 = cos(t40);
	t49 = t47 * t39;
	t38 = sin(t40);
	t54 = t45 * t38;
	t33 = t44 * t49 - t54;
	t31 = 0.1e1 / t33 ^ 2;
	t50 = t47 * t38;
	t53 = t45 * t39;
	t32 = t44 * t50 + t53;
	t58 = t31 * t32;
	t57 = t34 * t44;
	t43 = t46 ^ 2;
	t56 = 0.1e1 / t44 ^ 2 * t43;
	t36 = 0.1e1 / (t45 ^ 2 * t56 + 0.1e1);
	t55 = t45 * t36;
	t51 = t46 * t47;
	t48 = t32 ^ 2 * t31 + 0.1e1;
	t41 = 0.1e1 / t44;
	t30 = 0.1e1 / t33;
	t29 = (-0.1e1 - t56) * t55;
	t26 = 0.1e1 / t28;
	t25 = 0.1e1 / (t43 * t59 + 0.1e1);
	t24 = 0.1e1 / t48;
	t1 = [t41 * t36 * t51, 0, 0, t29, 0, 0; (t26 * t52 + (t35 * t41 * t43 * t55 + (-t36 + 0.1e1) * t46 * t34) * t46 * t59) * t25, 0, 0, (t44 * t26 + (-t45 * t57 + t35 * t46 + (t35 * t52 - t57) * t29) * t46 * t27) * t47 * t25, 0, 0; ((-t44 * t54 + t49) * t30 - (-t44 * t53 - t50) * t58) * t24, 0, 0, (t30 * t38 - t39 * t58) * t24 * t51, 0, t48 * t24;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end