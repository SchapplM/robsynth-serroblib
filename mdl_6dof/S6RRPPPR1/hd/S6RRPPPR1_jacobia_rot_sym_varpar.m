% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR1
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
%   Wie in S6RRPPPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPPR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (221->21), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->35)
	t38 = qJ(2) + pkin(9);
	t37 = cos(t38);
	t36 = sin(t38);
	t41 = sin(qJ(1));
	t47 = t41 * t36;
	t31 = atan2(-t47, -t37);
	t29 = sin(t31);
	t30 = cos(t31);
	t22 = -t29 * t47 - t30 * t37;
	t21 = 0.1e1 / t22 ^ 2;
	t42 = cos(qJ(1));
	t53 = t21 * t42 ^ 2;
	t40 = cos(pkin(10));
	t43 = t42 * t40;
	t39 = sin(pkin(10));
	t46 = t41 * t39;
	t28 = t37 * t43 + t46;
	t26 = 0.1e1 / t28 ^ 2;
	t44 = t42 * t39;
	t45 = t41 * t40;
	t27 = t37 * t44 - t45;
	t52 = t26 * t27;
	t51 = t29 * t37;
	t33 = t36 ^ 2;
	t50 = t33 / t37 ^ 2;
	t49 = t36 * t42;
	t32 = 0.1e1 / (t41 ^ 2 * t50 + 0.1e1);
	t48 = t41 * t32;
	t34 = 0.1e1 / t37;
	t25 = 0.1e1 / t28;
	t24 = 0.1e1 / (t27 ^ 2 * t26 + 0.1e1);
	t23 = (0.1e1 + t50) * t48;
	t20 = 0.1e1 / t22;
	t19 = 0.1e1 / (t33 * t53 + 0.1e1);
	t1 = [t34 * t32 * t49, t23, 0, 0, 0, 0; (-t20 * t47 - (-t30 * t33 * t34 * t48 + (t32 - 0.1e1) * t36 * t29) * t36 * t53) * t19, (t37 * t20 - (-t41 * t51 + t30 * t36 + (-t30 * t47 + t51) * t23) * t36 * t21) * t42 * t19, 0, 0, 0, 0; ((-t37 * t46 - t43) * t25 - (-t37 * t45 + t44) * t52) * t24, (-t25 * t39 + t40 * t52) * t24 * t49, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:09
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (227->21), mult. (329->51), div. (61->11), fcn. (494->9), ass. (0->36)
	t52 = qJ(2) + pkin(9);
	t49 = sin(t52);
	t69 = t49 ^ 2;
	t50 = cos(t52);
	t55 = cos(pkin(10));
	t57 = cos(qJ(1));
	t59 = t57 * t55;
	t54 = sin(pkin(10));
	t56 = sin(qJ(1));
	t62 = t56 * t54;
	t41 = t50 * t62 + t59;
	t63 = t49 * t54;
	t37 = atan2(-t41, t63);
	t34 = sin(t37);
	t35 = cos(t37);
	t33 = -t34 * t41 + t35 * t63;
	t32 = 0.1e1 / t33 ^ 2;
	t60 = t57 * t54;
	t61 = t56 * t55;
	t43 = t50 * t60 - t61;
	t68 = t32 * t43;
	t66 = t35 * t41;
	t65 = t43 ^ 2 * t32;
	t51 = 0.1e1 / t54;
	t64 = 0.1e1 / t49 * t51;
	t44 = t50 * t59 + t62;
	t40 = 0.1e1 / t44 ^ 2;
	t58 = t57 ^ 2 * t69 * t40;
	t48 = 0.1e1 / t69;
	t39 = 0.1e1 / t44;
	t38 = 0.1e1 / (0.1e1 + t58);
	t36 = 0.1e1 / (0.1e1 + t41 ^ 2 * t48 / t54 ^ 2);
	t31 = 0.1e1 / t33;
	t30 = (t41 * t48 * t50 * t51 + t56) * t36;
	t29 = 0.1e1 / (0.1e1 + t65);
	t1 = [-t43 * t36 * t64, t30, 0, 0, 0, 0; (-t41 * t31 - (-t34 + (t64 * t66 + t34) * t36) * t65) * t29, (t30 * t66 * t68 + (-t57 * t49 * t31 - (t35 * t50 + (-t30 + t56) * t34 * t49) * t68) * t54) * t29, 0, 0, 0, 0; (t56 * t39 + (-t50 * t61 + t60) * t57 * t40) * t49 * t38, (-t39 * t50 * t57 - t55 * t58) * t38, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:09
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (265->25), mult. (347->66), div. (52->9), fcn. (503->11), ass. (0->42)
	t58 = qJ(2) + pkin(9);
	t57 = cos(t58);
	t56 = sin(t58);
	t62 = sin(qJ(1));
	t70 = t62 * t56;
	t52 = atan2(t70, t57);
	t49 = sin(t52);
	t50 = cos(t52);
	t39 = t49 * t70 + t50 * t57;
	t38 = 0.1e1 / t39 ^ 2;
	t64 = cos(qJ(1));
	t76 = t38 * t64 ^ 2;
	t59 = sin(pkin(10));
	t67 = t64 * t59;
	t60 = cos(pkin(10));
	t68 = t62 * t60;
	t47 = t57 * t67 - t68;
	t66 = t64 * t60;
	t69 = t62 * t59;
	t48 = t57 * t66 + t69;
	t61 = sin(qJ(6));
	t63 = cos(qJ(6));
	t44 = t47 * t61 + t48 * t63;
	t42 = 0.1e1 / t44 ^ 2;
	t43 = -t47 * t63 + t48 * t61;
	t75 = t42 * t43;
	t74 = t49 * t57;
	t53 = t56 ^ 2;
	t73 = t53 / t57 ^ 2;
	t72 = t56 * t64;
	t51 = 0.1e1 / (t62 ^ 2 * t73 + 0.1e1);
	t71 = t62 * t51;
	t65 = t43 ^ 2 * t42 + 0.1e1;
	t54 = 0.1e1 / t57;
	t46 = -t57 * t68 + t67;
	t45 = -t57 * t69 - t66;
	t41 = 0.1e1 / t44;
	t40 = (0.1e1 + t73) * t71;
	t37 = 0.1e1 / t39;
	t36 = 0.1e1 / (t53 * t76 + 0.1e1);
	t35 = 0.1e1 / t65;
	t1 = [t54 * t51 * t72, t40, 0, 0, 0, 0; (t37 * t70 + (t50 * t53 * t54 * t71 + (-t51 + 0.1e1) * t56 * t49) * t56 * t76) * t36, (-t57 * t37 + (t62 * t74 - t50 * t56 + (t50 * t70 - t74) * t40) * t56 * t38) * t64 * t36, 0, 0, 0, 0; ((-t45 * t63 + t46 * t61) * t41 - (t45 * t61 + t46 * t63) * t75) * t35, ((t59 * t63 - t60 * t61) * t41 - (-t59 * t61 - t60 * t63) * t75) * t35 * t72, 0, 0, 0, t65 * t35;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end