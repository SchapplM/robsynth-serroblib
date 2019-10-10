% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR10
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
%   Wie in S6RPRPRR10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:02
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR10_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR10_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (76->19), mult. (197->54), div. (47->9), fcn. (297->9), ass. (0->35)
	t33 = sin(qJ(1));
	t47 = t33 ^ 2;
	t32 = sin(qJ(3));
	t34 = cos(qJ(3));
	t35 = cos(qJ(1));
	t36 = t35 * t34;
	t25 = atan2(-t36, t32);
	t23 = sin(t25);
	t24 = cos(t25);
	t16 = -t23 * t36 + t24 * t32;
	t15 = 0.1e1 / t16 ^ 2;
	t46 = t15 * t34;
	t30 = sin(pkin(10));
	t38 = t35 * t30;
	t31 = cos(pkin(10));
	t41 = t33 * t31;
	t22 = t32 * t41 + t38;
	t20 = 0.1e1 / t22 ^ 2;
	t37 = t35 * t31;
	t42 = t33 * t30;
	t21 = t32 * t42 - t37;
	t45 = t20 * t21;
	t44 = t23 * t32;
	t29 = t34 ^ 2;
	t43 = 0.1e1 / t32 ^ 2 * t29;
	t40 = t33 * t34;
	t26 = 0.1e1 / (t35 ^ 2 * t43 + 0.1e1);
	t39 = t35 * t26;
	t27 = 0.1e1 / t32;
	t19 = 0.1e1 / t22;
	t18 = (0.1e1 + t43) * t39;
	t17 = 0.1e1 / (t21 ^ 2 * t20 + 0.1e1);
	t14 = 0.1e1 / t16;
	t13 = 0.1e1 / (t47 * t29 * t15 + 0.1e1);
	t1 = [t27 * t26 * t40, 0, t18, 0, 0, 0; (-t14 * t36 + (-t24 * t27 * t29 * t39 + (-t26 + 0.1e1) * t34 * t23) * t47 * t46) * t13, 0, (t32 * t14 + (t35 * t44 + t24 * t34 + (-t24 * t36 - t44) * t18) * t46) * t33 * t13, 0, 0, 0; ((t32 * t38 + t41) * t19 - (t32 * t37 - t42) * t45) * t17, 0, (t19 * t30 - t31 * t45) * t17 * t40, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (134->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->37)
	t44 = sin(qJ(1));
	t59 = t44 ^ 2;
	t43 = sin(qJ(3));
	t45 = cos(qJ(3));
	t46 = cos(qJ(1));
	t48 = t46 * t45;
	t35 = atan2(-t48, t43);
	t33 = sin(t35);
	t34 = cos(t35);
	t27 = -t33 * t48 + t34 * t43;
	t26 = 0.1e1 / t27 ^ 2;
	t58 = t26 * t45;
	t39 = pkin(10) + qJ(5);
	t37 = sin(t39);
	t50 = t46 * t37;
	t38 = cos(t39);
	t53 = t44 * t38;
	t32 = t43 * t53 + t50;
	t30 = 0.1e1 / t32 ^ 2;
	t49 = t46 * t38;
	t54 = t44 * t37;
	t31 = t43 * t54 - t49;
	t57 = t30 * t31;
	t56 = t33 * t43;
	t42 = t45 ^ 2;
	t55 = 0.1e1 / t43 ^ 2 * t42;
	t52 = t44 * t45;
	t36 = 0.1e1 / (t46 ^ 2 * t55 + 0.1e1);
	t51 = t46 * t36;
	t47 = t31 ^ 2 * t30 + 0.1e1;
	t40 = 0.1e1 / t43;
	t29 = 0.1e1 / t32;
	t28 = (0.1e1 + t55) * t51;
	t25 = 0.1e1 / t27;
	t24 = 0.1e1 / t47;
	t23 = 0.1e1 / (t59 * t42 * t26 + 0.1e1);
	t1 = [t40 * t36 * t52, 0, t28, 0, 0, 0; (-t25 * t48 + (-t34 * t40 * t42 * t51 + (-t36 + 0.1e1) * t45 * t33) * t59 * t58) * t23, 0, (t43 * t25 + (t46 * t56 + t34 * t45 + (-t34 * t48 - t56) * t28) * t58) * t44 * t23, 0, 0, 0; ((t43 * t50 + t53) * t29 - (t43 * t49 - t54) * t57) * t24, 0, (t29 * t37 - t38 * t57) * t24 * t52, 0, t47 * t24, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (220->20), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->38)
	t57 = sin(qJ(1));
	t72 = t57 ^ 2;
	t56 = sin(qJ(3));
	t58 = cos(qJ(3));
	t59 = cos(qJ(1));
	t61 = t59 * t58;
	t48 = atan2(-t61, t56);
	t46 = sin(t48);
	t47 = cos(t48);
	t40 = -t46 * t61 + t47 * t56;
	t39 = 0.1e1 / t40 ^ 2;
	t71 = t39 * t58;
	t52 = pkin(10) + qJ(5) + qJ(6);
	t50 = sin(t52);
	t63 = t59 * t50;
	t51 = cos(t52);
	t66 = t57 * t51;
	t45 = t56 * t66 + t63;
	t43 = 0.1e1 / t45 ^ 2;
	t62 = t59 * t51;
	t67 = t57 * t50;
	t44 = t56 * t67 - t62;
	t70 = t43 * t44;
	t69 = t46 * t56;
	t55 = t58 ^ 2;
	t68 = 0.1e1 / t56 ^ 2 * t55;
	t65 = t57 * t58;
	t49 = 0.1e1 / (t59 ^ 2 * t68 + 0.1e1);
	t64 = t59 * t49;
	t60 = t44 ^ 2 * t43 + 0.1e1;
	t53 = 0.1e1 / t56;
	t42 = 0.1e1 / t45;
	t41 = (0.1e1 + t68) * t64;
	t38 = 0.1e1 / t40;
	t37 = 0.1e1 / (t72 * t55 * t39 + 0.1e1);
	t36 = 0.1e1 / t60;
	t35 = t60 * t36;
	t1 = [t53 * t49 * t65, 0, t41, 0, 0, 0; (-t38 * t61 + (-t47 * t53 * t55 * t64 + (-t49 + 0.1e1) * t58 * t46) * t72 * t71) * t37, 0, (t56 * t38 + (t59 * t69 + t47 * t58 + (-t47 * t61 - t69) * t41) * t71) * t57 * t37, 0, 0, 0; ((t56 * t63 + t66) * t42 - (t56 * t62 - t67) * t70) * t36, 0, (t42 * t50 - t51 * t70) * t36 * t65, 0, t35, t35;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end