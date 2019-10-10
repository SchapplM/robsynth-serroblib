% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR2
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
%   Wie in S6RRPPRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:37
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:31
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:31
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:31
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:31
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:31
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (221->21), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->35)
	t38 = qJ(2) + pkin(10);
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
	t40 = cos(pkin(11));
	t43 = t42 * t40;
	t39 = sin(pkin(11));
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
	% StartTime: 2019-10-10 09:37:31
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (286->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t52 = qJ(2) + pkin(10);
	t50 = cos(t52);
	t48 = sin(t52);
	t53 = sin(qJ(1));
	t59 = t53 * t48;
	t42 = atan2(-t59, -t50);
	t40 = sin(t42);
	t41 = cos(t42);
	t33 = -t40 * t59 - t41 * t50;
	t32 = 0.1e1 / t33 ^ 2;
	t54 = cos(qJ(1));
	t66 = t32 * t54 ^ 2;
	t51 = pkin(11) + qJ(5);
	t49 = cos(t51);
	t56 = t54 * t49;
	t47 = sin(t51);
	t60 = t53 * t47;
	t39 = t50 * t56 + t60;
	t37 = 0.1e1 / t39 ^ 2;
	t57 = t54 * t47;
	t58 = t53 * t49;
	t38 = t50 * t57 - t58;
	t65 = t37 * t38;
	t64 = t40 * t50;
	t44 = t48 ^ 2;
	t63 = t44 / t50 ^ 2;
	t62 = t48 * t54;
	t43 = 0.1e1 / (t53 ^ 2 * t63 + 0.1e1);
	t61 = t53 * t43;
	t55 = t38 ^ 2 * t37 + 0.1e1;
	t45 = 0.1e1 / t50;
	t36 = 0.1e1 / t39;
	t35 = (0.1e1 + t63) * t61;
	t34 = 0.1e1 / t55;
	t31 = 0.1e1 / t33;
	t30 = 0.1e1 / (t44 * t66 + 0.1e1);
	t1 = [t45 * t43 * t62, t35, 0, 0, 0, 0; (-t31 * t59 - (-t41 * t44 * t45 * t61 + (t43 - 0.1e1) * t48 * t40) * t48 * t66) * t30, (t50 * t31 - (-t53 * t64 + t41 * t48 + (-t41 * t59 + t64) * t35) * t48 * t32) * t54 * t30, 0, 0, 0, 0; ((-t50 * t60 - t56) * t36 - (-t50 * t58 + t57) * t65) * t34, (-t36 * t47 + t49 * t65) * t34 * t62, 0, 0, t55 * t34, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:31
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (379->22), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->38)
	t65 = qJ(2) + pkin(10);
	t63 = cos(t65);
	t62 = sin(t65);
	t66 = sin(qJ(1));
	t71 = t66 * t62;
	t55 = atan2(-t71, -t63);
	t53 = sin(t55);
	t54 = cos(t55);
	t47 = -t53 * t71 - t54 * t63;
	t46 = 0.1e1 / t47 ^ 2;
	t67 = cos(qJ(1));
	t79 = t46 * t67 ^ 2;
	t64 = pkin(11) + qJ(5) + qJ(6);
	t58 = cos(t64);
	t69 = t67 * t58;
	t57 = sin(t64);
	t73 = t66 * t57;
	t52 = t63 * t69 + t73;
	t50 = 0.1e1 / t52 ^ 2;
	t70 = t67 * t57;
	t72 = t66 * t58;
	t51 = t63 * t70 - t72;
	t78 = t50 * t51;
	t77 = t53 * t63;
	t59 = t62 ^ 2;
	t76 = t59 / t63 ^ 2;
	t75 = t62 * t67;
	t56 = 0.1e1 / (t66 ^ 2 * t76 + 0.1e1);
	t74 = t66 * t56;
	t68 = t51 ^ 2 * t50 + 0.1e1;
	t60 = 0.1e1 / t63;
	t49 = 0.1e1 / t52;
	t48 = (0.1e1 + t76) * t74;
	t45 = 0.1e1 / t47;
	t44 = 0.1e1 / t68;
	t43 = 0.1e1 / (t59 * t79 + 0.1e1);
	t42 = t68 * t44;
	t1 = [t60 * t56 * t75, t48, 0, 0, 0, 0; (-t45 * t71 - (-t54 * t59 * t60 * t74 + (t56 - 0.1e1) * t62 * t53) * t62 * t79) * t43, (t63 * t45 - (-t66 * t77 + t54 * t62 + (-t54 * t71 + t77) * t48) * t62 * t46) * t67 * t43, 0, 0, 0, 0; ((-t63 * t73 - t69) * t49 - (-t63 * t72 + t70) * t78) * t44, (-t49 * t57 + t58 * t78) * t44 * t75, 0, 0, t42, t42;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end