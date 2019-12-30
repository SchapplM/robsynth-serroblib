% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPP7
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
%   Wie in S5RRRPP7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:51
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRRPP7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_jacobia_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:51:06
	% EndTime: 2019-12-29 19:51:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:51:06
	% EndTime: 2019-12-29 19:51:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:51:06
	% EndTime: 2019-12-29 19:51:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:51:11
	% EndTime: 2019-12-29 19:51:11
	% DurationCPUTime: 0.15s
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
	% StartTime: 2019-12-29 19:51:06
	% EndTime: 2019-12-29 19:51:06
	% DurationCPUTime: 0.22s
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
	% StartTime: 2019-12-29 19:50:59
	% EndTime: 2019-12-29 19:50:59
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (64->19), mult. (227->55), div. (53->9), fcn. (337->9), ass. (0->33)
	t46 = cos(qJ(2));
	t43 = sin(qJ(2));
	t44 = sin(qJ(1));
	t51 = t44 * t43;
	t38 = atan2(t51, t46);
	t35 = sin(t38);
	t36 = cos(t38);
	t27 = t35 * t51 + t36 * t46;
	t26 = 0.1e1 / t27 ^ 2;
	t47 = cos(qJ(1));
	t57 = t26 * t47 ^ 2;
	t42 = sin(qJ(3));
	t45 = cos(qJ(3));
	t48 = t47 * t45;
	t34 = t44 * t42 + t46 * t48;
	t32 = 0.1e1 / t34 ^ 2;
	t49 = t47 * t42;
	t33 = t44 * t45 - t46 * t49;
	t56 = t33 ^ 2 * t32;
	t55 = t32 * t33;
	t39 = t43 ^ 2;
	t54 = t39 / t46 ^ 2;
	t53 = t43 * t47;
	t37 = 0.1e1 / (t44 ^ 2 * t54 + 0.1e1);
	t52 = t44 * t37;
	t50 = t44 * t46;
	t40 = 0.1e1 / t46;
	t31 = 0.1e1 / t34;
	t29 = (0.1e1 + t54) * t52;
	t28 = 0.1e1 / (0.1e1 + t56);
	t25 = 0.1e1 / t27;
	t24 = 0.1e1 / (t39 * t57 + 0.1e1);
	t1 = [t40 * t37 * t53, t29, 0, 0, 0; (t25 * t51 + (t36 * t39 * t40 * t52 + (-t37 + 0.1e1) * t43 * t35) * t43 * t57) * t24, (-t46 * t25 + (t35 * t50 - t36 * t43 + (-t35 * t46 + t36 * t51) * t29) * t43 * t26) * t47 * t24, 0, 0, 0; ((t42 * t50 + t48) * t31 - (-t45 * t50 + t49) * t55) * t28, (t31 * t42 + t45 * t55) * t28 * t53, (-t31 * t34 - t56) * t28, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end