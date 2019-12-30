% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRP8
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
%   Wie in S5RRRRP8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRRRP8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:38:47
	% EndTime: 2019-12-29 20:38:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:38:47
	% EndTime: 2019-12-29 20:38:47
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:38:47
	% EndTime: 2019-12-29 20:38:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:38:47
	% EndTime: 2019-12-29 20:38:47
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
	% StartTime: 2019-12-29 20:38:47
	% EndTime: 2019-12-29 20:38:47
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (181->21), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->35)
	t56 = cos(qJ(2));
	t54 = sin(qJ(2));
	t55 = sin(qJ(1));
	t62 = t55 * t54;
	t46 = atan2(-t62, -t56);
	t44 = sin(t46);
	t45 = cos(t46);
	t38 = -t44 * t62 - t45 * t56;
	t37 = 0.1e1 / t38 ^ 2;
	t57 = cos(qJ(1));
	t67 = t37 * t57 ^ 2;
	t53 = qJ(3) + qJ(4);
	t48 = sin(t53);
	t49 = cos(t53);
	t59 = t57 * t49;
	t43 = t55 * t48 + t56 * t59;
	t41 = 0.1e1 / t43 ^ 2;
	t60 = t57 * t48;
	t42 = -t55 * t49 + t56 * t60;
	t66 = t41 * t42;
	t50 = t54 ^ 2;
	t65 = t50 / t56 ^ 2;
	t64 = t54 * t57;
	t47 = 0.1e1 / (t55 ^ 2 * t65 + 0.1e1);
	t63 = t55 * t47;
	t61 = t55 * t56;
	t58 = t42 ^ 2 * t41 + 0.1e1;
	t51 = 0.1e1 / t56;
	t40 = 0.1e1 / t43;
	t39 = (0.1e1 + t65) * t63;
	t36 = 0.1e1 / t38;
	t35 = 0.1e1 / t58;
	t34 = 0.1e1 / (t50 * t67 + 0.1e1);
	t33 = t58 * t35;
	t1 = [t51 * t47 * t64, t39, 0, 0, 0; (-t36 * t62 - (-t45 * t50 * t51 * t63 + (t47 - 0.1e1) * t54 * t44) * t54 * t67) * t34, (t56 * t36 - (-t44 * t61 + t45 * t54 + (t44 * t56 - t45 * t62) * t39) * t54 * t37) * t57 * t34, 0, 0, 0; ((-t48 * t61 - t59) * t40 - (-t49 * t61 + t60) * t66) * t35, (-t40 * t48 + t49 * t66) * t35 * t64, t33, t33, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:38:47
	% EndTime: 2019-12-29 20:38:47
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (181->21), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->35)
	t57 = cos(qJ(2));
	t55 = sin(qJ(2));
	t56 = sin(qJ(1));
	t63 = t56 * t55;
	t47 = atan2(-t63, -t57);
	t45 = sin(t47);
	t46 = cos(t47);
	t39 = -t45 * t63 - t46 * t57;
	t38 = 0.1e1 / t39 ^ 2;
	t58 = cos(qJ(1));
	t68 = t38 * t58 ^ 2;
	t54 = qJ(3) + qJ(4);
	t49 = sin(t54);
	t50 = cos(t54);
	t60 = t58 * t50;
	t44 = t56 * t49 + t57 * t60;
	t42 = 0.1e1 / t44 ^ 2;
	t61 = t58 * t49;
	t43 = -t56 * t50 + t57 * t61;
	t67 = t42 * t43;
	t51 = t55 ^ 2;
	t66 = t51 / t57 ^ 2;
	t65 = t55 * t58;
	t48 = 0.1e1 / (t56 ^ 2 * t66 + 0.1e1);
	t64 = t56 * t48;
	t62 = t56 * t57;
	t59 = t43 ^ 2 * t42 + 0.1e1;
	t52 = 0.1e1 / t57;
	t41 = 0.1e1 / t44;
	t40 = (0.1e1 + t66) * t64;
	t37 = 0.1e1 / t39;
	t36 = 0.1e1 / t59;
	t35 = 0.1e1 / (t51 * t68 + 0.1e1);
	t34 = t59 * t36;
	t1 = [t52 * t48 * t65, t40, 0, 0, 0; (-t37 * t63 - (-t46 * t51 * t52 * t64 + (t48 - 0.1e1) * t55 * t45) * t55 * t68) * t35, (t57 * t37 - (-t45 * t62 + t46 * t55 + (t45 * t57 - t46 * t63) * t40) * t55 * t38) * t58 * t35, 0, 0, 0; ((-t49 * t62 - t60) * t41 - (-t50 * t62 + t61) * t67) * t36, (-t41 * t49 + t50 * t67) * t36 * t65, t34, t34, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end