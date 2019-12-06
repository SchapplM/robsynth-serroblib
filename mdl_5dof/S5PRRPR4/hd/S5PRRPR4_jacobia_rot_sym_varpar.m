% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRPR4
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
%   Wie in S5PRRPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRRPR4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:08
	% EndTime: 2019-12-05 16:24:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:08
	% EndTime: 2019-12-05 16:24:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:08
	% EndTime: 2019-12-05 16:24:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:08
	% EndTime: 2019-12-05 16:24:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (63->14), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->23)
	t31 = cos(qJ(2));
	t26 = sin(pkin(8));
	t29 = sin(qJ(2));
	t34 = t26 * t29;
	t22 = atan2(-t34, -t31);
	t20 = sin(t22);
	t36 = t20 * t31;
	t24 = t29 ^ 2;
	t35 = t24 / t31 ^ 2;
	t27 = cos(pkin(8));
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
	% StartTime: 2019-12-05 16:24:08
	% EndTime: 2019-12-05 16:24:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (93->15), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->24)
	t38 = cos(qJ(2));
	t35 = sin(pkin(8));
	t37 = sin(qJ(2));
	t41 = t35 * t37;
	t28 = atan2(-t41, -t38);
	t26 = sin(t28);
	t43 = t26 * t38;
	t33 = t37 ^ 2;
	t42 = t33 / t38 ^ 2;
	t36 = cos(pkin(8));
	t40 = t36 * t38;
	t32 = qJ(3) + pkin(9);
	t30 = sin(t32);
	t31 = cos(t32);
	t25 = t35 * t30 + t31 * t40;
	t23 = 0.1e1 / t25 ^ 2;
	t24 = t30 * t40 - t35 * t31;
	t39 = t24 ^ 2 * t23 + 0.1e1;
	t27 = cos(t28);
	t22 = (0.1e1 + t42) * t35 / (t35 ^ 2 * t42 + 0.1e1);
	t21 = -t26 * t41 - t27 * t38;
	t20 = 0.1e1 / t21 ^ 2;
	t19 = 0.1e1 / t39;
	t1 = [0, t22, 0, 0, 0; 0, (t38 / t21 - (-t35 * t43 + t27 * t37 + (-t27 * t41 + t43) * t22) * t37 * t20) * t36 / (t36 ^ 2 * t33 * t20 + 0.1e1), 0, 0, 0; 0, (-t30 / t25 + t31 * t24 * t23) * t37 * t36 * t19, t39 * t19, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:08
	% EndTime: 2019-12-05 16:24:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (161->15), mult. (162->34), div. (37->8), fcn. (227->9), ass. (0->25)
	t48 = cos(qJ(2));
	t45 = sin(pkin(8));
	t47 = sin(qJ(2));
	t51 = t45 * t47;
	t38 = atan2(-t51, -t48);
	t36 = sin(t38);
	t53 = t36 * t48;
	t43 = t47 ^ 2;
	t52 = t43 / t48 ^ 2;
	t46 = cos(pkin(8));
	t50 = t46 * t48;
	t42 = qJ(3) + pkin(9) + qJ(5);
	t40 = sin(t42);
	t41 = cos(t42);
	t35 = t45 * t40 + t41 * t50;
	t33 = 0.1e1 / t35 ^ 2;
	t34 = t40 * t50 - t45 * t41;
	t49 = t34 ^ 2 * t33 + 0.1e1;
	t37 = cos(t38);
	t32 = (0.1e1 + t52) * t45 / (t45 ^ 2 * t52 + 0.1e1);
	t31 = -t36 * t51 - t37 * t48;
	t30 = 0.1e1 / t31 ^ 2;
	t28 = 0.1e1 / t49;
	t27 = t49 * t28;
	t1 = [0, t32, 0, 0, 0; 0, (t48 / t31 - (-t45 * t53 + t37 * t47 + (-t37 * t51 + t53) * t32) * t47 * t30) * t46 / (t46 ^ 2 * t43 * t30 + 0.1e1), 0, 0, 0; 0, (-t40 / t35 + t41 * t34 * t33) * t47 * t46 * t28, t27, 0, t27;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end