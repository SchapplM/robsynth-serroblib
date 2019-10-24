% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRPP3
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
%   Wie in S5PRRPP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:29
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRRPP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:31
	% EndTime: 2019-10-24 10:29:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:31
	% EndTime: 2019-10-24 10:29:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:31
	% EndTime: 2019-10-24 10:29:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:31
	% EndTime: 2019-10-24 10:29:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (63->14), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->23)
	t31 = cos(qJ(2));
	t26 = sin(pkin(7));
	t29 = sin(qJ(2));
	t34 = t26 * t29;
	t22 = atan2(-t34, -t31);
	t20 = sin(t22);
	t36 = t20 * t31;
	t24 = t29 ^ 2;
	t35 = t24 / t31 ^ 2;
	t27 = cos(pkin(7));
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
	% StartTime: 2019-10-24 10:29:31
	% EndTime: 2019-10-24 10:29:32
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (139->25), mult. (431->69), div. (78->11), fcn. (644->11), ass. (0->39)
	t47 = sin(pkin(7));
	t49 = cos(pkin(7));
	t52 = cos(qJ(3));
	t50 = sin(qJ(3));
	t53 = cos(qJ(2));
	t57 = t50 * t53;
	t36 = t47 * t57 + t49 * t52;
	t51 = sin(qJ(2));
	t56 = t51 * t50;
	t35 = atan2(-t36, t56);
	t32 = sin(t35);
	t33 = cos(t35);
	t26 = -t32 * t36 + t33 * t56;
	t25 = 0.1e1 / t26 ^ 2;
	t39 = -t47 * t52 + t49 * t57;
	t62 = t25 * t39;
	t54 = t52 * t53;
	t40 = t47 * t50 + t49 * t54;
	t46 = sin(pkin(8));
	t48 = cos(pkin(8));
	t58 = t49 * t51;
	t31 = t40 * t48 + t46 * t58;
	t29 = 0.1e1 / t31 ^ 2;
	t30 = t40 * t46 - t48 * t58;
	t61 = t29 * t30;
	t59 = t33 * t36;
	t55 = t51 * t52;
	t45 = 0.1e1 / t51 ^ 2;
	t43 = 0.1e1 / t50 ^ 2;
	t42 = 0.1e1 / t50;
	t38 = t47 * t54 - t49 * t50;
	t34 = 0.1e1 / (t36 ^ 2 * t45 * t43 + 0.1e1);
	t28 = 0.1e1 / t31;
	t27 = 0.1e1 / (t30 ^ 2 * t29 + 0.1e1);
	t24 = 0.1e1 / t26;
	t23 = (t36 * t42 * t45 * t53 + t47) * t34;
	t22 = 0.1e1 / (t39 ^ 2 * t25 + 0.1e1);
	t21 = (t36 * t43 * t52 - t38 * t42) / t51 * t34;
	t1 = [0, t23, t21, 0, 0; 0, (t23 * t59 * t62 + (-t24 * t58 - (t33 * t53 + (-t23 + t47) * t51 * t32) * t62) * t50) * t22, (t40 * t24 - (t33 * t55 - t32 * t38 + (-t32 * t56 - t59) * t21) * t62) * t22, 0, 0; 0, ((-t46 * t55 - t48 * t53) * t28 - (t46 * t53 - t48 * t55) * t61) * t27 * t49, (-t28 * t46 + t48 * t61) * t39 * t27, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:31
	% EndTime: 2019-10-24 10:29:32
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (256->26), mult. (752->70), div. (60->9), fcn. (1075->11), ass. (0->42)
	t62 = sin(pkin(8));
	t63 = sin(pkin(7));
	t65 = cos(pkin(7));
	t66 = sin(qJ(3));
	t68 = cos(qJ(3));
	t69 = cos(qJ(2));
	t71 = t68 * t69;
	t64 = cos(pkin(8));
	t67 = sin(qJ(2));
	t73 = t67 * t64;
	t48 = (t63 * t71 - t65 * t66) * t62 - t63 * t73;
	t72 = t67 * t68;
	t57 = t62 * t72 + t64 * t69;
	t45 = atan2(-t48, t57);
	t41 = sin(t45);
	t42 = cos(t45);
	t40 = -t41 * t48 + t42 * t57;
	t39 = 0.1e1 / t40 ^ 2;
	t59 = t63 * t66 + t65 * t71;
	t50 = t59 * t62 - t65 * t73;
	t78 = t39 * t50;
	t54 = 0.1e1 / t57 ^ 2;
	t77 = t48 * t54;
	t51 = t62 * t65 * t67 + t59 * t64;
	t47 = 0.1e1 / t51 ^ 2;
	t74 = t66 * t69;
	t58 = t63 * t68 - t65 * t74;
	t76 = t58 ^ 2 * t47;
	t75 = t66 * t67;
	t70 = -t41 * t57 - t42 * t48;
	t60 = t62 * t71 - t73;
	t56 = -t63 * t74 - t65 * t68;
	t53 = 0.1e1 / t57;
	t52 = t57 * t63;
	t46 = 0.1e1 / t51;
	t44 = 0.1e1 / (t48 ^ 2 * t54 + 0.1e1);
	t43 = 0.1e1 / (0.1e1 + t76);
	t38 = 0.1e1 / t40;
	t37 = 0.1e1 / (t39 * t50 ^ 2 + 0.1e1);
	t36 = (-t53 * t56 - t75 * t77) * t62 * t44;
	t35 = (t52 * t53 + t60 * t77) * t44;
	t1 = [0, t35, t36, 0, 0; 0, (-(t35 * t70 + t41 * t52 + t42 * t60) * t78 - t57 * t38 * t65) * t37, (t58 * t62 * t38 - ((-t41 * t56 - t42 * t75) * t62 + t70 * t36) * t78) * t37, 0, 0; 0, (t46 * t75 - (t62 * t69 - t64 * t72) * t58 * t47) * t43 * t65, (-t46 * t59 - t64 * t76) * t43, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end