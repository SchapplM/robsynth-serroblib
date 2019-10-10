% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR6
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
%   Wie in S6RPRRPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:32
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:12
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (238->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
	t43 = pkin(10) + qJ(3);
	t42 = cos(t43);
	t41 = sin(t43);
	t45 = sin(qJ(1));
	t53 = t45 * t41;
	t36 = atan2(-t53, -t42);
	t32 = sin(t36);
	t33 = cos(t36);
	t27 = -t32 * t53 - t33 * t42;
	t26 = 0.1e1 / t27 ^ 2;
	t47 = cos(qJ(1));
	t59 = t26 * t47 ^ 2;
	t46 = cos(qJ(4));
	t49 = t47 * t46;
	t44 = sin(qJ(4));
	t52 = t45 * t44;
	t35 = t42 * t49 + t52;
	t31 = 0.1e1 / t35 ^ 2;
	t50 = t47 * t44;
	t51 = t45 * t46;
	t34 = t42 * t50 - t51;
	t58 = t31 * t34;
	t57 = t32 * t42;
	t38 = t41 ^ 2;
	t56 = t38 / t42 ^ 2;
	t55 = t41 * t47;
	t37 = 0.1e1 / (t45 ^ 2 * t56 + 0.1e1);
	t54 = t45 * t37;
	t48 = t34 ^ 2 * t31 + 0.1e1;
	t39 = 0.1e1 / t42;
	t30 = 0.1e1 / t35;
	t29 = 0.1e1 / t48;
	t28 = (0.1e1 + t56) * t54;
	t25 = 0.1e1 / t27;
	t24 = 0.1e1 / (t38 * t59 + 0.1e1);
	t1 = [t39 * t37 * t55, 0, t28, 0, 0, 0; (-t25 * t53 - (-t33 * t38 * t39 * t54 + (t37 - 0.1e1) * t41 * t32) * t41 * t59) * t24, 0, (t42 * t25 - (-t45 * t57 + t33 * t41 + (-t33 * t53 + t57) * t28) * t41 * t26) * t47 * t24, 0, 0, 0; ((-t42 * t52 - t49) * t30 - (-t42 * t51 + t50) * t58) * t29, 0, (-t30 * t44 + t46 * t58) * t29 * t55, t48 * t29, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:12
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (286->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t52 = pkin(10) + qJ(3);
	t50 = cos(t52);
	t48 = sin(t52);
	t54 = sin(qJ(1));
	t61 = t54 * t48;
	t43 = atan2(-t61, -t50);
	t41 = sin(t43);
	t42 = cos(t43);
	t34 = -t41 * t61 - t42 * t50;
	t33 = 0.1e1 / t34 ^ 2;
	t55 = cos(qJ(1));
	t67 = t33 * t55 ^ 2;
	t53 = qJ(4) + pkin(11);
	t51 = cos(t53);
	t57 = t55 * t51;
	t49 = sin(t53);
	t60 = t54 * t49;
	t40 = t50 * t57 + t60;
	t38 = 0.1e1 / t40 ^ 2;
	t58 = t55 * t49;
	t59 = t54 * t51;
	t39 = t50 * t58 - t59;
	t66 = t38 * t39;
	t65 = t41 * t50;
	t45 = t48 ^ 2;
	t64 = t45 / t50 ^ 2;
	t63 = t48 * t55;
	t44 = 0.1e1 / (t54 ^ 2 * t64 + 0.1e1);
	t62 = t54 * t44;
	t56 = t39 ^ 2 * t38 + 0.1e1;
	t46 = 0.1e1 / t50;
	t37 = 0.1e1 / t40;
	t36 = (0.1e1 + t64) * t62;
	t35 = 0.1e1 / t56;
	t32 = 0.1e1 / t34;
	t31 = 0.1e1 / (t45 * t67 + 0.1e1);
	t1 = [t46 * t44 * t63, 0, t36, 0, 0, 0; (-t32 * t61 - (-t42 * t45 * t46 * t62 + (t44 - 0.1e1) * t48 * t41) * t48 * t67) * t31, 0, (t50 * t32 - (-t54 * t65 + t42 * t48 + (-t42 * t61 + t65) * t36) * t48 * t33) * t55 * t31, 0, 0, 0; ((-t50 * t60 - t57) * t37 - (-t50 * t59 + t58) * t66) * t35, 0, (-t37 * t49 + t51 * t66) * t35 * t63, t56 * t35, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:12
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
	t61 = qJ(4) + pkin(11) + qJ(6);
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
	t1 = [t57 * t53 * t72, 0, t45, 0, 0, 0; (-t42 * t68 - (-t51 * t56 * t57 * t71 + (t53 - 0.1e1) * t59 * t50) * t59 * t76) * t40, 0, (t60 * t42 - (-t63 * t74 + t51 * t59 + (-t51 * t68 + t74) * t45) * t59 * t43) * t64 * t40, 0, 0, 0; ((-t60 * t70 - t66) * t46 - (-t60 * t69 + t67) * t75) * t41, 0, (-t46 * t54 + t55 * t75) * t41 * t72, t39, 0, t39;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end