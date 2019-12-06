% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRP6
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
%   Wie in S5PRRRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRRRP6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:53:03
	% EndTime: 2019-12-05 16:53:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:53:03
	% EndTime: 2019-12-05 16:53:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:53:03
	% EndTime: 2019-12-05 16:53:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:53:03
	% EndTime: 2019-12-05 16:53:03
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-12-05 16:53:03
	% EndTime: 2019-12-05 16:53:03
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (117->15), mult. (162->34), div. (37->8), fcn. (227->9), ass. (0->25)
	t46 = cos(qJ(2));
	t43 = sin(pkin(8));
	t45 = sin(qJ(2));
	t49 = t43 * t45;
	t36 = atan2(-t49, -t46);
	t34 = sin(t36);
	t51 = t34 * t46;
	t40 = t45 ^ 2;
	t50 = t40 / t46 ^ 2;
	t44 = cos(pkin(8));
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
	% StartTime: 2019-12-05 16:53:03
	% EndTime: 2019-12-05 16:53:03
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (479->23), mult. (539->54), div. (114->11), fcn. (816->9), ass. (0->37)
	t70 = sin(qJ(2));
	t80 = t70 ^ 2;
	t67 = qJ(3) + qJ(4);
	t61 = sin(t67);
	t62 = cos(t67);
	t69 = cos(pkin(8));
	t68 = sin(pkin(8));
	t71 = cos(qJ(2));
	t76 = t68 * t71;
	t54 = t61 * t76 + t69 * t62;
	t73 = t70 * t61;
	t51 = atan2(-t54, t73);
	t48 = sin(t51);
	t49 = cos(t51);
	t46 = -t48 * t54 + t49 * t73;
	t45 = 0.1e1 / t46 ^ 2;
	t74 = t69 * t71;
	t57 = t61 * t74 - t68 * t62;
	t79 = t45 * t57;
	t77 = t49 * t54;
	t75 = t69 * t70;
	t58 = t68 * t61 + t62 * t74;
	t53 = 0.1e1 / t58 ^ 2;
	t72 = t69 ^ 2 * t80 * t53;
	t66 = 0.1e1 / t80;
	t60 = 0.1e1 / t61 ^ 2;
	t59 = 0.1e1 / t61;
	t56 = -t69 * t61 + t62 * t76;
	t52 = 0.1e1 / (0.1e1 + t72);
	t50 = 0.1e1 / (t54 ^ 2 * t66 * t60 + 0.1e1);
	t47 = t57 * t53 * t52 * t75;
	t44 = 0.1e1 / t46;
	t43 = (t54 * t59 * t66 * t71 + t68) * t50;
	t42 = 0.1e1 / (t57 ^ 2 * t45 + 0.1e1);
	t41 = (t54 * t60 * t62 - t56 * t59) / t70 * t50;
	t40 = (t58 * t44 - (t49 * t70 * t62 - t48 * t56 + (-t48 * t73 - t77) * t41) * t79) * t42;
	t1 = [0, t43, t41, t41, 0; 0, (t43 * t77 * t79 + (-t44 * t75 - (t49 * t71 + (-t43 + t68) * t70 * t48) * t79) * t61) * t42, t40, t40, 0; 0, (-0.1e1 / t58 * t74 - t62 * t72) * t52, -t47, -t47, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end