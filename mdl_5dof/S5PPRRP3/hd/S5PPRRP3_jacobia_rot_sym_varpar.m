% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRRP3
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
%   Wie in S5PPRRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:20
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PPRRP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:03
	% EndTime: 2019-10-24 10:20:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:04
	% EndTime: 2019-10-24 10:20:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:03
	% EndTime: 2019-10-24 10:20:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:03
	% EndTime: 2019-10-24 10:20:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 1, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:04
	% EndTime: 2019-10-24 10:20:04
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (89->17), mult. (268->41), div. (47->11), fcn. (395->11), ass. (0->31)
	t42 = sin(pkin(8));
	t45 = cos(pkin(7));
	t56 = t42 * t45;
	t47 = sin(qJ(3));
	t55 = t42 * t47;
	t43 = sin(pkin(7));
	t54 = t43 * t47;
	t49 = cos(qJ(3));
	t53 = t43 * t49;
	t52 = t45 * t47;
	t51 = t45 * t49;
	t44 = cos(pkin(8));
	t38 = t44 * t51 + t54;
	t46 = sin(qJ(4));
	t48 = cos(qJ(4));
	t29 = t38 * t48 + t46 * t56;
	t27 = 0.1e1 / t29 ^ 2;
	t28 = t38 * t46 - t48 * t56;
	t50 = t28 ^ 2 * t27 + 0.1e1;
	t41 = 0.1e1 / t47 ^ 2;
	t37 = t44 * t52 - t53;
	t36 = t44 * t53 - t52;
	t34 = t44 * t54 + t51;
	t33 = atan2(-t34, t55);
	t31 = cos(t33);
	t30 = sin(t33);
	t26 = 0.1e1 / t50;
	t25 = -t30 * t34 + t31 * t55;
	t24 = 0.1e1 / t25 ^ 2;
	t22 = (-t36 / t47 + t49 * t34 * t41) / t42 / (0.1e1 + t34 ^ 2 / t42 ^ 2 * t41);
	t1 = [0, 0, t22, 0, 0; 0, 0, (t38 / t25 - (t31 * t42 * t49 - t30 * t36 + (-t30 * t55 - t31 * t34) * t22) * t37 * t24) / (t37 ^ 2 * t24 + 0.1e1), 0, 0; 0, 0, (-t46 / t29 + t48 * t28 * t27) * t37 * t26, t50 * t26, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:04
	% EndTime: 2019-10-24 10:20:04
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (256->25), mult. (739->64), div. (57->9), fcn. (1061->11), ass. (0->45)
	t61 = cos(pkin(8));
	t62 = cos(pkin(7));
	t64 = sin(qJ(3));
	t69 = t62 * t64;
	t60 = sin(pkin(7));
	t66 = cos(qJ(3));
	t70 = t60 * t66;
	t53 = t61 * t70 - t69;
	t63 = sin(qJ(4));
	t59 = sin(pkin(8));
	t65 = cos(qJ(4));
	t73 = t59 * t65;
	t44 = t53 * t63 - t60 * t73;
	t72 = t59 * t66;
	t56 = t61 * t65 + t63 * t72;
	t42 = atan2(-t44, t56);
	t38 = sin(t42);
	t39 = cos(t42);
	t37 = -t38 * t44 + t39 * t56;
	t36 = 0.1e1 / t37 ^ 2;
	t68 = t62 * t66;
	t71 = t60 * t64;
	t55 = t61 * t68 + t71;
	t47 = t55 * t63 - t62 * t73;
	t78 = t36 * t47;
	t51 = 0.1e1 / t56 ^ 2;
	t77 = t44 * t51;
	t75 = t59 * t63;
	t48 = t55 * t65 + t62 * t75;
	t43 = 0.1e1 / t48 ^ 2;
	t54 = -t61 * t69 + t70;
	t76 = t54 ^ 2 * t43;
	t74 = t59 * t64;
	t67 = -t38 * t56 - t39 * t44;
	t57 = -t61 * t63 + t65 * t72;
	t52 = -t61 * t71 - t68;
	t50 = 0.1e1 / t56;
	t46 = t53 * t65 + t60 * t75;
	t41 = 0.1e1 / (t44 ^ 2 * t51 + 0.1e1);
	t40 = 0.1e1 / (0.1e1 + t76);
	t35 = 0.1e1 / t37;
	t34 = 0.1e1 / (t47 ^ 2 * t36 + 0.1e1);
	t33 = (-t50 * t52 - t74 * t77) * t63 * t41;
	t32 = (-t46 * t50 + t57 * t77) * t41;
	t1 = [0, 0, t33, t32, 0; 0, 0, (t54 * t63 * t35 - ((-t38 * t52 - t39 * t74) * t63 + t67 * t33) * t78) * t34, (t48 * t35 - (t67 * t32 - t38 * t46 + t39 * t57) * t78) * t34, 0; 0, 0, (-t55 / t48 - t65 * t76) * t40, t47 * t54 * t43 * t40, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end