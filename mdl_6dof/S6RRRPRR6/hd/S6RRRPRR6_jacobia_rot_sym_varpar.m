% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR6
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
%   Wie in S6RRRPRR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:02
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.10s
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
	t1 = [t34 * t32 * t48, t24, 0, 0, 0, 0; (-t20 * t46 - (-t30 * t33 * t34 * t47 + (t32 - 0.1e1) * t37 * t29) * t37 * t51) * t19, (t40 * t20 - (-t29 * t45 + t30 * t37 + (t29 * t40 - t30 * t46) * t24) * t37 * t21) * t41 * t19, 0, 0, 0, 0; ((-t36 * t45 - t43) * t25 - (-t39 * t45 + t44) * t50) * t23, (-t25 * t36 + t39 * t50) * t23 * t48, t42 * t23, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (157->21), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->34)
	t48 = cos(qJ(2));
	t46 = sin(qJ(2));
	t47 = sin(qJ(1));
	t54 = t47 * t46;
	t38 = atan2(-t54, -t48);
	t36 = sin(t38);
	t37 = cos(t38);
	t30 = -t36 * t54 - t37 * t48;
	t29 = 0.1e1 / t30 ^ 2;
	t49 = cos(qJ(1));
	t59 = t29 * t49 ^ 2;
	t42 = qJ(3) + pkin(11);
	t40 = sin(t42);
	t41 = cos(t42);
	t51 = t49 * t41;
	t35 = t47 * t40 + t48 * t51;
	t33 = 0.1e1 / t35 ^ 2;
	t52 = t49 * t40;
	t34 = -t47 * t41 + t48 * t52;
	t58 = t33 * t34;
	t43 = t46 ^ 2;
	t57 = t43 / t48 ^ 2;
	t56 = t46 * t49;
	t39 = 0.1e1 / (t47 ^ 2 * t57 + 0.1e1);
	t55 = t47 * t39;
	t53 = t47 * t48;
	t50 = t34 ^ 2 * t33 + 0.1e1;
	t44 = 0.1e1 / t48;
	t32 = 0.1e1 / t35;
	t31 = (0.1e1 + t57) * t55;
	t28 = 0.1e1 / t30;
	t27 = 0.1e1 / t50;
	t26 = 0.1e1 / (t43 * t59 + 0.1e1);
	t1 = [t44 * t39 * t56, t31, 0, 0, 0, 0; (-t28 * t54 - (-t37 * t43 * t44 * t55 + (t39 - 0.1e1) * t46 * t36) * t46 * t59) * t26, (t48 * t28 - (-t36 * t53 + t37 * t46 + (t36 * t48 - t37 * t54) * t31) * t46 * t29) * t49 * t26, 0, 0, 0, 0; ((-t40 * t53 - t51) * t32 - (-t41 * t53 + t52) * t58) * t27, (-t32 * t40 + t41 * t58) * t27 * t56, t50 * t27, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (243->21), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->35)
	t59 = cos(qJ(2));
	t57 = sin(qJ(2));
	t58 = sin(qJ(1));
	t65 = t58 * t57;
	t49 = atan2(-t65, -t59);
	t47 = sin(t49);
	t48 = cos(t49);
	t41 = -t47 * t65 - t48 * t59;
	t40 = 0.1e1 / t41 ^ 2;
	t60 = cos(qJ(1));
	t70 = t40 * t60 ^ 2;
	t53 = qJ(3) + pkin(11) + qJ(5);
	t51 = sin(t53);
	t52 = cos(t53);
	t62 = t60 * t52;
	t46 = t58 * t51 + t59 * t62;
	t44 = 0.1e1 / t46 ^ 2;
	t63 = t60 * t51;
	t45 = -t58 * t52 + t59 * t63;
	t69 = t44 * t45;
	t54 = t57 ^ 2;
	t68 = t54 / t59 ^ 2;
	t67 = t57 * t60;
	t50 = 0.1e1 / (t58 ^ 2 * t68 + 0.1e1);
	t66 = t58 * t50;
	t64 = t58 * t59;
	t61 = t45 ^ 2 * t44 + 0.1e1;
	t55 = 0.1e1 / t59;
	t43 = 0.1e1 / t46;
	t42 = (0.1e1 + t68) * t66;
	t39 = 0.1e1 / t41;
	t38 = 0.1e1 / (t54 * t70 + 0.1e1);
	t37 = 0.1e1 / t61;
	t36 = t61 * t37;
	t1 = [t55 * t50 * t67, t42, 0, 0, 0, 0; (-t39 * t65 - (-t48 * t54 * t55 * t66 + (t50 - 0.1e1) * t57 * t47) * t57 * t70) * t38, (t59 * t39 - (-t47 * t64 + t48 * t57 + (t47 * t59 - t48 * t65) * t42) * t57 * t40) * t60 * t38, 0, 0, 0, 0; ((-t51 * t64 - t62) * t43 - (-t52 * t64 + t63) * t69) * t37, (-t43 * t51 + t52 * t69) * t37 * t67, t36, 0, t36, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (357->21), mult. (278->55), div. (62->9), fcn. (402->9), ass. (0->35)
	t65 = cos(qJ(2));
	t63 = sin(qJ(2));
	t64 = sin(qJ(1));
	t71 = t64 * t63;
	t55 = atan2(-t71, -t65);
	t53 = sin(t55);
	t54 = cos(t55);
	t47 = -t53 * t71 - t54 * t65;
	t46 = 0.1e1 / t47 ^ 2;
	t66 = cos(qJ(1));
	t76 = t46 * t66 ^ 2;
	t59 = qJ(3) + pkin(11) + qJ(5) + qJ(6);
	t57 = sin(t59);
	t58 = cos(t59);
	t68 = t66 * t58;
	t52 = t64 * t57 + t65 * t68;
	t50 = 0.1e1 / t52 ^ 2;
	t69 = t66 * t57;
	t51 = -t64 * t58 + t65 * t69;
	t75 = t50 * t51;
	t60 = t63 ^ 2;
	t74 = t60 / t65 ^ 2;
	t73 = t63 * t66;
	t56 = 0.1e1 / (t64 ^ 2 * t74 + 0.1e1);
	t72 = t64 * t56;
	t70 = t64 * t65;
	t67 = t51 ^ 2 * t50 + 0.1e1;
	t61 = 0.1e1 / t65;
	t49 = 0.1e1 / t52;
	t48 = (0.1e1 + t74) * t72;
	t45 = 0.1e1 / t47;
	t44 = 0.1e1 / (t60 * t76 + 0.1e1);
	t43 = 0.1e1 / t67;
	t42 = t67 * t43;
	t1 = [t61 * t56 * t73, t48, 0, 0, 0, 0; (-t45 * t71 - (-t54 * t60 * t61 * t72 + (t56 - 0.1e1) * t63 * t53) * t63 * t76) * t44, (t65 * t45 - (-t53 * t70 + t54 * t63 + (t53 * t65 - t54 * t71) * t48) * t63 * t46) * t66 * t44, 0, 0, 0, 0; ((-t57 * t70 - t68) * t49 - (-t58 * t70 + t69) * t75) * t43, (-t49 * t57 + t58 * t75) * t43 * t73, t42, 0, t42, t42;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end