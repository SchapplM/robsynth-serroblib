% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR2
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
%   Wie in S6RRPPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:19
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPPR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (192->17), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->28)
	t36 = cos(qJ(1));
	t34 = t36 ^ 2;
	t32 = qJ(2) + pkin(9);
	t31 = cos(t32);
	t30 = sin(t32);
	t35 = sin(qJ(1));
	t39 = t35 * t30;
	t24 = atan2(-t39, -t31);
	t22 = sin(t24);
	t23 = cos(t24);
	t20 = -t22 * t39 - t23 * t31;
	t19 = 0.1e1 / t20 ^ 2;
	t45 = t19 * t30;
	t44 = t22 * t31;
	t27 = t30 ^ 2;
	t37 = t31 ^ 2;
	t43 = t27 / t37;
	t42 = t30 * t36;
	t38 = t35 ^ 2;
	t41 = 0.1e1 / t38 * t34;
	t25 = 0.1e1 / (t38 * t43 + 0.1e1);
	t40 = t35 * t25;
	t28 = 0.1e1 / t31;
	t26 = 0.1e1 / (t37 * t41 + 0.1e1);
	t21 = (0.1e1 + t43) * t40;
	t18 = 0.1e1 / t20;
	t17 = 0.1e1 / (t34 * t27 * t19 + 0.1e1);
	t1 = [t28 * t25 * t42, t21, 0, 0, 0, 0; (-t18 * t39 - (-t23 * t27 * t28 * t40 + (t25 - 0.1e1) * t30 * t22) * t34 * t45) * t17, (t31 * t18 - (-t35 * t44 + t23 * t30 + (-t23 * t39 + t44) * t21) * t45) * t36 * t17, 0, 0, 0, 0; (-0.1e1 - t41) * t31 * t26, -0.1e1 / t35 * t26 * t42, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:53
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (199->21), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->35)
	t38 = qJ(2) + pkin(9);
	t36 = sin(t38);
	t37 = cos(t38);
	t41 = sin(qJ(1));
	t47 = t41 * t37;
	t31 = atan2(-t47, t36);
	t29 = sin(t31);
	t30 = cos(t31);
	t22 = -t29 * t47 + t30 * t36;
	t21 = 0.1e1 / t22 ^ 2;
	t42 = cos(qJ(1));
	t53 = t21 * t42 ^ 2;
	t39 = sin(pkin(10));
	t44 = t42 * t39;
	t40 = cos(pkin(10));
	t45 = t41 * t40;
	t28 = t36 * t44 + t45;
	t26 = 0.1e1 / t28 ^ 2;
	t43 = t42 * t40;
	t46 = t41 * t39;
	t27 = -t36 * t43 + t46;
	t52 = t26 * t27;
	t51 = t29 * t36;
	t35 = t37 ^ 2;
	t50 = 0.1e1 / t36 ^ 2 * t35;
	t49 = t37 * t42;
	t32 = 0.1e1 / (t41 ^ 2 * t50 + 0.1e1);
	t48 = t41 * t32;
	t33 = 0.1e1 / t36;
	t25 = 0.1e1 / t28;
	t24 = 0.1e1 / (t27 ^ 2 * t26 + 0.1e1);
	t23 = (0.1e1 + t50) * t48;
	t20 = 0.1e1 / t22;
	t19 = 0.1e1 / (t35 * t53 + 0.1e1);
	t1 = [-t33 * t32 * t49, t23, 0, 0, 0, 0; (-t20 * t47 - (t30 * t33 * t35 * t48 + (t32 - 0.1e1) * t37 * t29) * t37 * t53) * t19, (-t36 * t20 - (t41 * t51 + t30 * t37 + (-t30 * t47 - t51) * t23) * t37 * t21) * t42 * t19, 0, 0, 0, 0; ((t36 * t45 + t44) * t25 - (-t36 * t46 + t43) * t52) * t24, (-t25 * t40 - t39 * t52) * t24 * t49, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:19:53
	% EndTime: 2019-10-10 09:19:54
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (264->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t55 = qJ(2) + pkin(9);
	t51 = sin(t55);
	t53 = cos(t55);
	t56 = sin(qJ(1));
	t61 = t56 * t53;
	t45 = atan2(-t61, t51);
	t43 = sin(t45);
	t44 = cos(t45);
	t36 = -t43 * t61 + t44 * t51;
	t35 = 0.1e1 / t36 ^ 2;
	t57 = cos(qJ(1));
	t69 = t35 * t57 ^ 2;
	t54 = pkin(10) + qJ(6);
	t50 = sin(t54);
	t60 = t57 * t50;
	t52 = cos(t54);
	t62 = t56 * t52;
	t42 = t51 * t60 + t62;
	t40 = 0.1e1 / t42 ^ 2;
	t59 = t57 * t52;
	t63 = t56 * t50;
	t41 = -t51 * t59 + t63;
	t68 = t40 * t41;
	t67 = t43 * t51;
	t49 = t53 ^ 2;
	t66 = 0.1e1 / t51 ^ 2 * t49;
	t65 = t53 * t57;
	t46 = 0.1e1 / (t56 ^ 2 * t66 + 0.1e1);
	t64 = t56 * t46;
	t58 = t41 ^ 2 * t40 + 0.1e1;
	t47 = 0.1e1 / t51;
	t39 = 0.1e1 / t42;
	t38 = (0.1e1 + t66) * t64;
	t37 = 0.1e1 / t58;
	t34 = 0.1e1 / t36;
	t33 = 0.1e1 / (t49 * t69 + 0.1e1);
	t1 = [-t47 * t46 * t65, t38, 0, 0, 0, 0; (-t34 * t61 - (t44 * t47 * t49 * t64 + (t46 - 0.1e1) * t53 * t43) * t53 * t69) * t33, (-t51 * t34 - (t56 * t67 + t44 * t53 + (-t44 * t61 - t67) * t38) * t53 * t35) * t57 * t33, 0, 0, 0, 0; ((t51 * t62 + t60) * t39 - (-t51 * t63 + t59) * t68) * t37, (-t39 * t52 - t50 * t68) * t37 * t65, 0, 0, 0, t58 * t37;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end