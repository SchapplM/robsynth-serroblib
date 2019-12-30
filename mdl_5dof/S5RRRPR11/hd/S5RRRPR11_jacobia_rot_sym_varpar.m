% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPR11
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
%   Wie in S5RRRPR11_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRRPR11_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR11_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:25
	% EndTime: 2019-12-29 20:16:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:13
	% EndTime: 2019-12-29 20:16:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:19
	% EndTime: 2019-12-29 20:16:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:18
	% EndTime: 2019-12-29 20:16:19
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (109->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->33)
	t40 = cos(qJ(2));
	t37 = sin(qJ(2));
	t38 = sin(qJ(1));
	t46 = t38 * t37;
	t31 = atan2(-t46, -t40);
	t29 = sin(t31);
	t30 = cos(t31);
	t22 = -t29 * t46 - t30 * t40;
	t21 = 0.1e1 / t22 ^ 2;
	t41 = cos(qJ(1));
	t51 = t21 * t41 ^ 2;
	t36 = sin(qJ(3));
	t39 = cos(qJ(3));
	t43 = t41 * t39;
	t28 = t38 * t36 + t40 * t43;
	t26 = 0.1e1 / t28 ^ 2;
	t44 = t41 * t36;
	t27 = -t38 * t39 + t40 * t44;
	t50 = t26 * t27;
	t33 = t37 ^ 2;
	t49 = t33 / t40 ^ 2;
	t48 = t37 * t41;
	t32 = 0.1e1 / (t38 ^ 2 * t49 + 0.1e1);
	t47 = t38 * t32;
	t45 = t38 * t40;
	t42 = t27 ^ 2 * t26 + 0.1e1;
	t34 = 0.1e1 / t40;
	t25 = 0.1e1 / t28;
	t24 = (0.1e1 + t49) * t47;
	t23 = 0.1e1 / t42;
	t20 = 0.1e1 / t22;
	t19 = 0.1e1 / (t33 * t51 + 0.1e1);
	t1 = [t34 * t32 * t48, t24, 0, 0, 0; (-t20 * t46 - (-t30 * t33 * t34 * t47 + (t32 - 0.1e1) * t37 * t29) * t37 * t51) * t19, (t40 * t20 - (-t29 * t45 + t30 * t37 + (t29 * t40 - t30 * t46) * t24) * t37 * t21) * t41 * t19, 0, 0, 0; ((-t36 * t45 - t43) * t25 - (-t39 * t45 + t44) * t50) * t23, (-t25 * t36 + t39 * t50) * t23 * t48, t42 * t23, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:25
	% EndTime: 2019-12-29 20:16:26
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (159->26), mult. (509->70), div. (100->11), fcn. (770->9), ass. (0->39)
	t52 = sin(qJ(2));
	t68 = t52 ^ 2;
	t51 = sin(qJ(3));
	t54 = cos(qJ(3));
	t56 = cos(qJ(1));
	t58 = t56 * t54;
	t53 = sin(qJ(1));
	t55 = cos(qJ(2));
	t60 = t53 * t55;
	t38 = t51 * t60 + t58;
	t62 = t52 * t51;
	t35 = atan2(-t38, t62);
	t31 = sin(t35);
	t32 = cos(t35);
	t30 = -t31 * t38 + t32 * t62;
	t29 = 0.1e1 / t30 ^ 2;
	t59 = t56 * t51;
	t41 = -t53 * t54 + t55 * t59;
	t67 = t29 * t41;
	t65 = t32 * t38;
	t64 = t41 ^ 2 * t29;
	t45 = 0.1e1 / t51;
	t48 = 0.1e1 / t52;
	t63 = t45 * t48;
	t61 = t52 * t56;
	t42 = t53 * t51 + t55 * t58;
	t37 = 0.1e1 / t42 ^ 2;
	t57 = t56 ^ 2 * t68 * t37;
	t49 = 0.1e1 / t68;
	t46 = 0.1e1 / t51 ^ 2;
	t40 = t54 * t60 - t59;
	t36 = 0.1e1 / t42;
	t34 = 0.1e1 / (t38 ^ 2 * t49 * t46 + 0.1e1);
	t33 = 0.1e1 / (0.1e1 + t57);
	t28 = 0.1e1 / t30;
	t27 = (t38 * t45 * t49 * t55 + t53) * t34;
	t26 = 0.1e1 / (0.1e1 + t64);
	t25 = (t38 * t46 * t54 - t40 * t45) * t48 * t34;
	t1 = [-t41 * t34 * t63, t27, t25, 0, 0; (-t38 * t28 - (-t31 + (t63 * t65 + t31) * t34) * t64) * t26, (t27 * t65 * t67 + (-t28 * t61 - (t32 * t55 + (-t27 + t53) * t52 * t31) * t67) * t51) * t26, (t42 * t28 - (t32 * t52 * t54 - t31 * t40 + (-t31 * t62 - t65) * t25) * t67) * t26, 0, 0; (-t37 * t40 * t56 + t36 * t53) * t52 * t33, (-t36 * t55 * t56 - t54 * t57) * t33, -t41 * t37 * t33 * t61, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:25
	% EndTime: 2019-12-29 20:16:26
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (141->25), mult. (425->68), div. (58->9), fcn. (611->11), ass. (0->40)
	t60 = sin(qJ(1));
	t62 = cos(qJ(3));
	t63 = cos(qJ(2));
	t58 = sin(qJ(3));
	t64 = cos(qJ(1));
	t67 = t64 * t58;
	t47 = -t60 * t62 + t63 * t67;
	t66 = t64 * t62;
	t48 = t60 * t58 + t63 * t66;
	t57 = sin(qJ(5));
	t61 = cos(qJ(5));
	t40 = t47 * t57 + t48 * t61;
	t38 = 0.1e1 / t40 ^ 2;
	t39 = -t47 * t61 + t48 * t57;
	t75 = t38 * t39;
	t74 = t39 ^ 2 * t38;
	t59 = sin(qJ(2));
	t69 = t60 * t59;
	t52 = atan2(t69, t63);
	t49 = sin(t52);
	t50 = cos(t52);
	t43 = t49 * t69 + t50 * t63;
	t42 = 0.1e1 / t43 ^ 2;
	t73 = t42 * t64 ^ 2;
	t54 = t59 ^ 2;
	t72 = t54 / t63 ^ 2;
	t71 = t59 * t64;
	t51 = 0.1e1 / (t60 ^ 2 * t72 + 0.1e1);
	t70 = t60 * t51;
	t68 = t60 * t63;
	t65 = 0.1e1 + t74;
	t55 = 0.1e1 / t63;
	t46 = -t62 * t68 + t67;
	t45 = -t58 * t68 - t66;
	t44 = (0.1e1 + t72) * t70;
	t41 = 0.1e1 / t43;
	t37 = 0.1e1 / t40;
	t36 = 0.1e1 / (t54 * t73 + 0.1e1);
	t35 = 0.1e1 / t65;
	t1 = [t55 * t51 * t71, t44, 0, 0, 0; (t41 * t69 + (t50 * t54 * t55 * t70 + (-t51 + 0.1e1) * t59 * t49) * t59 * t73) * t36, (-t63 * t41 + (t49 * t68 - t50 * t59 + (-t49 * t63 + t50 * t69) * t44) * t59 * t42) * t64 * t36, 0, 0, 0; ((-t45 * t61 + t46 * t57) * t37 - (t45 * t57 + t46 * t61) * t75) * t35, ((-t57 * t62 + t58 * t61) * t37 - (-t57 * t58 - t61 * t62) * t75) * t35 * t71, (-t37 * t40 - t74) * t35, 0, t65 * t35;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end