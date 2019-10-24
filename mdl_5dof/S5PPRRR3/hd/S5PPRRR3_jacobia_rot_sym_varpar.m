% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRRR3
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
%   Wie in S5PPRRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:21
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PPRRR3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:17
	% EndTime: 2019-10-24 10:21:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:17
	% EndTime: 2019-10-24 10:21:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:17
	% EndTime: 2019-10-24 10:21:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:17
	% EndTime: 2019-10-24 10:21:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 1, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:17
	% EndTime: 2019-10-24 10:21:18
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (89->17), mult. (268->41), div. (47->11), fcn. (395->11), ass. (0->31)
	t42 = sin(pkin(9));
	t45 = cos(pkin(8));
	t56 = t42 * t45;
	t47 = sin(qJ(3));
	t55 = t42 * t47;
	t43 = sin(pkin(8));
	t54 = t43 * t47;
	t49 = cos(qJ(3));
	t53 = t43 * t49;
	t52 = t45 * t47;
	t51 = t45 * t49;
	t44 = cos(pkin(9));
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
	% StartTime: 2019-10-24 10:21:18
	% EndTime: 2019-10-24 10:21:18
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (150->18), mult. (316->41), div. (52->11), fcn. (458->11), ass. (0->33)
	t61 = sin(pkin(9));
	t64 = cos(pkin(8));
	t73 = t61 * t64;
	t65 = sin(qJ(3));
	t72 = t61 * t65;
	t62 = sin(pkin(8));
	t71 = t62 * t65;
	t66 = cos(qJ(3));
	t70 = t62 * t66;
	t69 = t64 * t65;
	t68 = t64 * t66;
	t63 = cos(pkin(9));
	t54 = t63 * t68 + t71;
	t60 = qJ(4) + qJ(5);
	t56 = sin(t60);
	t57 = cos(t60);
	t45 = t54 * t57 + t56 * t73;
	t43 = 0.1e1 / t45 ^ 2;
	t44 = t54 * t56 - t57 * t73;
	t67 = t44 ^ 2 * t43 + 0.1e1;
	t59 = 0.1e1 / t65 ^ 2;
	t53 = t63 * t69 - t70;
	t52 = t63 * t70 - t69;
	t50 = t63 * t71 + t68;
	t49 = atan2(-t50, t72);
	t47 = cos(t49);
	t46 = sin(t49);
	t42 = 0.1e1 / t67;
	t41 = -t46 * t50 + t47 * t72;
	t40 = 0.1e1 / t41 ^ 2;
	t38 = (-t52 / t65 + t66 * t50 * t59) / t61 / (0.1e1 + t50 ^ 2 / t61 ^ 2 * t59);
	t37 = t67 * t42;
	t1 = [0, 0, t38, 0, 0; 0, 0, (t54 / t41 - (t47 * t61 * t66 - t46 * t52 + (-t46 * t72 - t47 * t50) * t38) * t53 * t40) / (t53 ^ 2 * t40 + 0.1e1), 0, 0; 0, 0, (-t56 / t45 + t57 * t44 * t43) * t53 * t42, t37, t37;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end