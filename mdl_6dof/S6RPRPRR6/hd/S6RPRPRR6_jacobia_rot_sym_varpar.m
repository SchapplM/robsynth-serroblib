% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR6
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
%   Wie in S6RPRPRR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:54:55
	% EndTime: 2019-10-10 00:54:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:54:55
	% EndTime: 2019-10-10 00:54:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:54:55
	% EndTime: 2019-10-10 00:54:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:54:55
	% EndTime: 2019-10-10 00:54:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:54:55
	% EndTime: 2019-10-10 00:54:55
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (221->21), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->35)
	t35 = pkin(10) + qJ(3);
	t34 = cos(t35);
	t33 = sin(t35);
	t38 = sin(qJ(1));
	t44 = t38 * t33;
	t28 = atan2(-t44, -t34);
	t26 = sin(t28);
	t27 = cos(t28);
	t19 = -t26 * t44 - t27 * t34;
	t18 = 0.1e1 / t19 ^ 2;
	t39 = cos(qJ(1));
	t50 = t18 * t39 ^ 2;
	t37 = cos(pkin(11));
	t40 = t39 * t37;
	t36 = sin(pkin(11));
	t43 = t38 * t36;
	t25 = t34 * t40 + t43;
	t23 = 0.1e1 / t25 ^ 2;
	t41 = t39 * t36;
	t42 = t38 * t37;
	t24 = t34 * t41 - t42;
	t49 = t23 * t24;
	t48 = t26 * t34;
	t30 = t33 ^ 2;
	t47 = t30 / t34 ^ 2;
	t46 = t33 * t39;
	t29 = 0.1e1 / (t38 ^ 2 * t47 + 0.1e1);
	t45 = t38 * t29;
	t31 = 0.1e1 / t34;
	t22 = 0.1e1 / t25;
	t21 = 0.1e1 / (t24 ^ 2 * t23 + 0.1e1);
	t20 = (0.1e1 + t47) * t45;
	t17 = 0.1e1 / t19;
	t16 = 0.1e1 / (t30 * t50 + 0.1e1);
	t1 = [t31 * t29 * t46, 0, t20, 0, 0, 0; (-t17 * t44 - (-t27 * t30 * t31 * t45 + (t29 - 0.1e1) * t33 * t26) * t33 * t50) * t16, 0, (t34 * t17 - (-t38 * t48 + t27 * t33 + (-t27 * t44 + t48) * t20) * t33 * t18) * t39 * t16, 0, 0, 0; ((-t34 * t43 - t40) * t22 - (-t34 * t42 + t41) * t49) * t21, 0, (-t22 * t36 + t37 * t49) * t21 * t46, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:54:55
	% EndTime: 2019-10-10 00:54:55
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (286->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t49 = pkin(10) + qJ(3);
	t47 = cos(t49);
	t45 = sin(t49);
	t50 = sin(qJ(1));
	t56 = t50 * t45;
	t39 = atan2(-t56, -t47);
	t37 = sin(t39);
	t38 = cos(t39);
	t30 = -t37 * t56 - t38 * t47;
	t29 = 0.1e1 / t30 ^ 2;
	t51 = cos(qJ(1));
	t63 = t29 * t51 ^ 2;
	t48 = pkin(11) + qJ(5);
	t46 = cos(t48);
	t53 = t51 * t46;
	t44 = sin(t48);
	t57 = t50 * t44;
	t36 = t47 * t53 + t57;
	t34 = 0.1e1 / t36 ^ 2;
	t54 = t51 * t44;
	t55 = t50 * t46;
	t35 = t47 * t54 - t55;
	t62 = t34 * t35;
	t61 = t37 * t47;
	t41 = t45 ^ 2;
	t60 = t41 / t47 ^ 2;
	t59 = t45 * t51;
	t40 = 0.1e1 / (t50 ^ 2 * t60 + 0.1e1);
	t58 = t50 * t40;
	t52 = t35 ^ 2 * t34 + 0.1e1;
	t42 = 0.1e1 / t47;
	t33 = 0.1e1 / t36;
	t32 = (0.1e1 + t60) * t58;
	t31 = 0.1e1 / t52;
	t28 = 0.1e1 / t30;
	t27 = 0.1e1 / (t41 * t63 + 0.1e1);
	t1 = [t42 * t40 * t59, 0, t32, 0, 0, 0; (-t28 * t56 - (-t38 * t41 * t42 * t58 + (t40 - 0.1e1) * t45 * t37) * t45 * t63) * t27, 0, (t47 * t28 - (-t50 * t61 + t38 * t45 + (-t38 * t56 + t61) * t32) * t45 * t29) * t51 * t27, 0, 0, 0; ((-t47 * t57 - t53) * t33 - (-t47 * t55 + t54) * t62) * t31, 0, (-t33 * t44 + t46 * t62) * t31 * t59, 0, t52 * t31, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:54:55
	% EndTime: 2019-10-10 00:54:55
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (379->22), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->38)
	t62 = pkin(10) + qJ(3);
	t60 = cos(t62);
	t59 = sin(t62);
	t63 = sin(qJ(1));
	t68 = t63 * t59;
	t52 = atan2(-t68, -t60);
	t50 = sin(t52);
	t51 = cos(t52);
	t44 = -t50 * t68 - t51 * t60;
	t43 = 0.1e1 / t44 ^ 2;
	t64 = cos(qJ(1));
	t76 = t43 * t64 ^ 2;
	t61 = pkin(11) + qJ(5) + qJ(6);
	t55 = cos(t61);
	t66 = t64 * t55;
	t54 = sin(t61);
	t70 = t63 * t54;
	t49 = t60 * t66 + t70;
	t47 = 0.1e1 / t49 ^ 2;
	t67 = t64 * t54;
	t69 = t63 * t55;
	t48 = t60 * t67 - t69;
	t75 = t47 * t48;
	t74 = t50 * t60;
	t56 = t59 ^ 2;
	t73 = t56 / t60 ^ 2;
	t72 = t59 * t64;
	t53 = 0.1e1 / (t63 ^ 2 * t73 + 0.1e1);
	t71 = t63 * t53;
	t65 = t48 ^ 2 * t47 + 0.1e1;
	t57 = 0.1e1 / t60;
	t46 = 0.1e1 / t49;
	t45 = (0.1e1 + t73) * t71;
	t42 = 0.1e1 / t44;
	t41 = 0.1e1 / t65;
	t40 = 0.1e1 / (t56 * t76 + 0.1e1);
	t39 = t65 * t41;
	t1 = [t57 * t53 * t72, 0, t45, 0, 0, 0; (-t42 * t68 - (-t51 * t56 * t57 * t71 + (t53 - 0.1e1) * t59 * t50) * t59 * t76) * t40, 0, (t60 * t42 - (-t63 * t74 + t51 * t59 + (-t51 * t68 + t74) * t45) * t59 * t43) * t64 * t40, 0, 0, 0; ((-t60 * t70 - t66) * t46 - (-t60 * t69 + t67) * t75) * t41, 0, (-t46 * t54 + t55 * t75) * t41 * t72, 0, t39, t39;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end