% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR5
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
%   Wie in S6RPRPPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:22
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (192->17), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->28)
	t33 = cos(qJ(1));
	t31 = t33 ^ 2;
	t29 = pkin(9) + qJ(3);
	t28 = cos(t29);
	t27 = sin(t29);
	t32 = sin(qJ(1));
	t36 = t32 * t27;
	t21 = atan2(-t36, -t28);
	t19 = sin(t21);
	t20 = cos(t21);
	t17 = -t19 * t36 - t20 * t28;
	t16 = 0.1e1 / t17 ^ 2;
	t42 = t16 * t27;
	t41 = t19 * t28;
	t24 = t27 ^ 2;
	t34 = t28 ^ 2;
	t40 = t24 / t34;
	t39 = t27 * t33;
	t35 = t32 ^ 2;
	t38 = 0.1e1 / t35 * t31;
	t22 = 0.1e1 / (t35 * t40 + 0.1e1);
	t37 = t32 * t22;
	t25 = 0.1e1 / t28;
	t23 = 0.1e1 / (t34 * t38 + 0.1e1);
	t18 = (0.1e1 + t40) * t37;
	t15 = 0.1e1 / t17;
	t14 = 0.1e1 / (t31 * t24 * t16 + 0.1e1);
	t1 = [t25 * t22 * t39, 0, t18, 0, 0, 0; (-t15 * t36 - (-t20 * t24 * t25 * t37 + (t22 - 0.1e1) * t27 * t19) * t31 * t42) * t14, 0, (t28 * t15 - (-t32 * t41 + t20 * t27 + (-t20 * t36 + t41) * t18) * t42) * t33 * t14, 0, 0, 0; (-0.1e1 - t38) * t28 * t23, 0, -0.1e1 / t32 * t23 * t39, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (199->21), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->35)
	t35 = pkin(9) + qJ(3);
	t33 = sin(t35);
	t34 = cos(t35);
	t38 = sin(qJ(1));
	t44 = t38 * t34;
	t28 = atan2(-t44, t33);
	t26 = sin(t28);
	t27 = cos(t28);
	t19 = -t26 * t44 + t27 * t33;
	t18 = 0.1e1 / t19 ^ 2;
	t39 = cos(qJ(1));
	t50 = t18 * t39 ^ 2;
	t36 = sin(pkin(10));
	t41 = t39 * t36;
	t37 = cos(pkin(10));
	t42 = t38 * t37;
	t25 = t33 * t41 + t42;
	t23 = 0.1e1 / t25 ^ 2;
	t40 = t39 * t37;
	t43 = t38 * t36;
	t24 = -t33 * t40 + t43;
	t49 = t23 * t24;
	t48 = t26 * t33;
	t32 = t34 ^ 2;
	t47 = 0.1e1 / t33 ^ 2 * t32;
	t46 = t34 * t39;
	t29 = 0.1e1 / (t38 ^ 2 * t47 + 0.1e1);
	t45 = t38 * t29;
	t30 = 0.1e1 / t33;
	t22 = 0.1e1 / t25;
	t21 = 0.1e1 / (t24 ^ 2 * t23 + 0.1e1);
	t20 = (0.1e1 + t47) * t45;
	t17 = 0.1e1 / t19;
	t16 = 0.1e1 / (t32 * t50 + 0.1e1);
	t1 = [-t30 * t29 * t46, 0, t20, 0, 0, 0; (-t17 * t44 - (t27 * t30 * t32 * t45 + (t29 - 0.1e1) * t34 * t26) * t34 * t50) * t16, 0, (-t33 * t17 - (t38 * t48 + t27 * t34 + (-t27 * t44 - t48) * t20) * t34 * t18) * t39 * t16, 0, 0, 0; ((t33 * t42 + t41) * t22 - (-t33 * t43 + t40) * t49) * t21, 0, (-t22 * t37 - t36 * t49) * t21 * t46, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:18
	% EndTime: 2019-10-10 00:22:18
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (264->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t52 = pkin(9) + qJ(3);
	t48 = sin(t52);
	t50 = cos(t52);
	t53 = sin(qJ(1));
	t58 = t53 * t50;
	t42 = atan2(-t58, t48);
	t40 = sin(t42);
	t41 = cos(t42);
	t33 = -t40 * t58 + t41 * t48;
	t32 = 0.1e1 / t33 ^ 2;
	t54 = cos(qJ(1));
	t66 = t32 * t54 ^ 2;
	t51 = pkin(10) + qJ(6);
	t47 = sin(t51);
	t57 = t54 * t47;
	t49 = cos(t51);
	t59 = t53 * t49;
	t39 = t48 * t57 + t59;
	t37 = 0.1e1 / t39 ^ 2;
	t56 = t54 * t49;
	t60 = t53 * t47;
	t38 = -t48 * t56 + t60;
	t65 = t37 * t38;
	t64 = t40 * t48;
	t46 = t50 ^ 2;
	t63 = 0.1e1 / t48 ^ 2 * t46;
	t62 = t50 * t54;
	t43 = 0.1e1 / (t53 ^ 2 * t63 + 0.1e1);
	t61 = t53 * t43;
	t55 = t38 ^ 2 * t37 + 0.1e1;
	t44 = 0.1e1 / t48;
	t36 = 0.1e1 / t39;
	t35 = (0.1e1 + t63) * t61;
	t34 = 0.1e1 / t55;
	t31 = 0.1e1 / t33;
	t30 = 0.1e1 / (t46 * t66 + 0.1e1);
	t1 = [-t44 * t43 * t62, 0, t35, 0, 0, 0; (-t31 * t58 - (t41 * t44 * t46 * t61 + (t43 - 0.1e1) * t50 * t40) * t50 * t66) * t30, 0, (-t48 * t31 - (t53 * t64 + t41 * t50 + (-t41 * t58 - t64) * t35) * t50 * t32) * t54 * t30, 0, 0, 0; ((t48 * t59 + t57) * t36 - (-t48 * t60 + t56) * t65) * t34, 0, (-t36 * t49 - t47 * t65) * t34 * t62, 0, 0, t55 * t34;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end