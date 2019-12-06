% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRR7
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
%   Wie in S5PRRRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRRRR7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:13:42
	% EndTime: 2019-12-05 17:13:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:13:42
	% EndTime: 2019-12-05 17:13:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:13:42
	% EndTime: 2019-12-05 17:13:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:13:41
	% EndTime: 2019-12-05 17:13:42
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (63->14), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->23)
	t31 = cos(qJ(2));
	t26 = sin(pkin(9));
	t29 = sin(qJ(2));
	t34 = t26 * t29;
	t22 = atan2(-t34, -t31);
	t20 = sin(t22);
	t36 = t20 * t31;
	t24 = t29 ^ 2;
	t35 = t24 / t31 ^ 2;
	t27 = cos(pkin(9));
	t33 = t27 * t31;
	t28 = sin(qJ(3));
	t30 = cos(qJ(3));
	t19 = t26 * t28 + t30 * t33;
	t17 = 0.1e1 / t19 ^ 2;
	t18 = -t26 * t30 + t28 * t33;
	t32 = t18 ^ 2 * t17 + 0.1e1;
	t21 = cos(t22);
	t16 = (0.1e1 + t35) * t26 / (t26 ^ 2 * t35 + 0.1e1);
	t15 = 0.1e1 / t32;
	t14 = -t20 * t34 - t21 * t31;
	t13 = 0.1e1 / t14 ^ 2;
	t1 = [0, t16, 0, 0, 0; 0, (t31 / t14 - (-t26 * t36 + t21 * t29 + (-t21 * t34 + t36) * t16) * t29 * t13) * t27 / (t27 ^ 2 * t24 * t13 + 0.1e1), 0, 0, 0; 0, (-t28 / t19 + t30 * t18 * t17) * t29 * t27 * t15, t32 * t15, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:13:42
	% EndTime: 2019-12-05 17:13:42
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (117->15), mult. (162->34), div. (37->8), fcn. (227->9), ass. (0->25)
	t46 = cos(qJ(2));
	t43 = sin(pkin(9));
	t45 = sin(qJ(2));
	t49 = t43 * t45;
	t36 = atan2(-t49, -t46);
	t34 = sin(t36);
	t51 = t34 * t46;
	t40 = t45 ^ 2;
	t50 = t40 / t46 ^ 2;
	t44 = cos(pkin(9));
	t48 = t44 * t46;
	t42 = qJ(3) + qJ(4);
	t38 = sin(t42);
	t39 = cos(t42);
	t33 = t43 * t38 + t39 * t48;
	t31 = 0.1e1 / t33 ^ 2;
	t32 = t38 * t48 - t43 * t39;
	t47 = t32 ^ 2 * t31 + 0.1e1;
	t35 = cos(t36);
	t30 = (0.1e1 + t50) * t43 / (t43 ^ 2 * t50 + 0.1e1);
	t29 = -t34 * t49 - t35 * t46;
	t28 = 0.1e1 / t29 ^ 2;
	t27 = 0.1e1 / t47;
	t25 = t47 * t27;
	t1 = [0, t30, 0, 0, 0; 0, (t46 / t29 - (-t43 * t51 + t35 * t45 + (-t35 * t49 + t51) * t30) * t45 * t28) * t44 / (t44 ^ 2 * t40 * t28 + 0.1e1), 0, 0, 0; 0, (-t38 / t33 + t39 * t32 * t31) * t45 * t44 * t27, t25, t25, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:13:42
	% EndTime: 2019-12-05 17:13:42
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (199->15), mult. (189->34), div. (42->8), fcn. (262->9), ass. (0->25)
	t51 = cos(qJ(2));
	t48 = sin(pkin(9));
	t50 = sin(qJ(2));
	t54 = t48 * t50;
	t41 = atan2(-t54, -t51);
	t39 = sin(t41);
	t56 = t39 * t51;
	t46 = t50 ^ 2;
	t55 = t46 / t51 ^ 2;
	t49 = cos(pkin(9));
	t53 = t49 * t51;
	t45 = qJ(3) + qJ(4) + qJ(5);
	t43 = sin(t45);
	t44 = cos(t45);
	t38 = t48 * t43 + t44 * t53;
	t36 = 0.1e1 / t38 ^ 2;
	t37 = t43 * t53 - t48 * t44;
	t52 = t37 ^ 2 * t36 + 0.1e1;
	t40 = cos(t41);
	t35 = (0.1e1 + t55) * t48 / (t48 ^ 2 * t55 + 0.1e1);
	t34 = -t39 * t54 - t40 * t51;
	t33 = 0.1e1 / t34 ^ 2;
	t31 = 0.1e1 / t52;
	t30 = t52 * t31;
	t1 = [0, t35, 0, 0, 0; 0, (t51 / t34 - (-t48 * t56 + t40 * t50 + (-t40 * t54 + t56) * t35) * t50 * t33) * t49 / (t49 ^ 2 * t46 * t33 + 0.1e1), 0, 0, 0; 0, (-t43 / t38 + t44 * t37 * t36) * t50 * t49 * t31, t30, t30, t30;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end