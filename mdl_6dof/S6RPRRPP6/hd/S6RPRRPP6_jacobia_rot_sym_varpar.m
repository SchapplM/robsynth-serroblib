% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP6
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
%   Wie in S6RPRRPP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:18
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:18:23
	% EndTime: 2019-10-10 01:18:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:18:23
	% EndTime: 2019-10-10 01:18:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:18:23
	% EndTime: 2019-10-10 01:18:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:18:23
	% EndTime: 2019-10-10 01:18:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:18:23
	% EndTime: 2019-10-10 01:18:24
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (86->19), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->36)
	t40 = sin(qJ(1));
	t56 = t40 ^ 2;
	t39 = sin(qJ(3));
	t42 = cos(qJ(3));
	t43 = cos(qJ(1));
	t45 = t43 * t42;
	t33 = atan2(-t45, t39);
	t31 = sin(t33);
	t32 = cos(t33);
	t24 = -t31 * t45 + t32 * t39;
	t23 = 0.1e1 / t24 ^ 2;
	t55 = t23 * t42;
	t38 = sin(qJ(4));
	t47 = t43 * t38;
	t41 = cos(qJ(4));
	t50 = t40 * t41;
	t30 = t39 * t50 + t47;
	t28 = 0.1e1 / t30 ^ 2;
	t46 = t43 * t41;
	t51 = t40 * t38;
	t29 = t39 * t51 - t46;
	t54 = t28 * t29;
	t53 = t31 * t39;
	t37 = t42 ^ 2;
	t52 = 0.1e1 / t39 ^ 2 * t37;
	t49 = t40 * t42;
	t34 = 0.1e1 / (t43 ^ 2 * t52 + 0.1e1);
	t48 = t43 * t34;
	t44 = t29 ^ 2 * t28 + 0.1e1;
	t35 = 0.1e1 / t39;
	t27 = 0.1e1 / t30;
	t26 = (0.1e1 + t52) * t48;
	t25 = 0.1e1 / t44;
	t22 = 0.1e1 / t24;
	t21 = 0.1e1 / (t56 * t37 * t23 + 0.1e1);
	t1 = [t35 * t34 * t49, 0, t26, 0, 0, 0; (-t22 * t45 + (-t32 * t35 * t37 * t48 + (-t34 + 0.1e1) * t42 * t31) * t56 * t55) * t21, 0, (t39 * t22 + (t43 * t53 + t32 * t42 + (-t32 * t45 - t53) * t26) * t55) * t40 * t21, 0, 0, 0; ((t39 * t47 + t50) * t27 - (t39 * t46 - t51) * t54) * t25, 0, (t27 * t38 - t41 * t54) * t25 * t49, t44 * t25, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:18:23
	% EndTime: 2019-10-10 01:18:24
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (134->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->37)
	t46 = sin(qJ(1));
	t61 = t46 ^ 2;
	t45 = sin(qJ(3));
	t47 = cos(qJ(3));
	t48 = cos(qJ(1));
	t50 = t48 * t47;
	t37 = atan2(-t50, t45);
	t35 = sin(t37);
	t36 = cos(t37);
	t29 = -t35 * t50 + t36 * t45;
	t28 = 0.1e1 / t29 ^ 2;
	t60 = t28 * t47;
	t41 = qJ(4) + pkin(9);
	t39 = sin(t41);
	t52 = t48 * t39;
	t40 = cos(t41);
	t55 = t46 * t40;
	t34 = t45 * t55 + t52;
	t32 = 0.1e1 / t34 ^ 2;
	t51 = t48 * t40;
	t56 = t46 * t39;
	t33 = t45 * t56 - t51;
	t59 = t32 * t33;
	t58 = t35 * t45;
	t44 = t47 ^ 2;
	t57 = 0.1e1 / t45 ^ 2 * t44;
	t54 = t46 * t47;
	t38 = 0.1e1 / (t48 ^ 2 * t57 + 0.1e1);
	t53 = t48 * t38;
	t49 = t33 ^ 2 * t32 + 0.1e1;
	t42 = 0.1e1 / t45;
	t31 = 0.1e1 / t34;
	t30 = (0.1e1 + t57) * t53;
	t27 = 0.1e1 / t29;
	t26 = 0.1e1 / t49;
	t25 = 0.1e1 / (t61 * t44 * t28 + 0.1e1);
	t1 = [t42 * t38 * t54, 0, t30, 0, 0, 0; (-t27 * t50 + (-t36 * t42 * t44 * t53 + (-t38 + 0.1e1) * t47 * t35) * t61 * t60) * t25, 0, (t45 * t27 + (t48 * t58 + t36 * t47 + (-t36 * t50 - t58) * t30) * t60) * t46 * t25, 0, 0, 0; ((t45 * t52 + t55) * t31 - (t45 * t51 - t56) * t59) * t26, 0, (t31 * t39 - t40 * t59) * t26 * t54, t49 * t26, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:18:23
	% EndTime: 2019-10-10 01:18:24
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (447->25), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->41)
	t59 = cos(qJ(3));
	t73 = t59 ^ 2;
	t57 = sin(qJ(3));
	t52 = qJ(4) + pkin(9);
	t50 = sin(t52);
	t60 = cos(qJ(1));
	t63 = t60 * t50;
	t51 = cos(t52);
	t58 = sin(qJ(1));
	t66 = t58 * t51;
	t45 = t57 * t63 + t66;
	t64 = t59 * t50;
	t39 = atan2(t45, t64);
	t36 = sin(t39);
	t37 = cos(t39);
	t35 = t36 * t45 + t37 * t64;
	t34 = 0.1e1 / t35 ^ 2;
	t62 = t60 * t51;
	t67 = t58 * t50;
	t43 = t57 * t67 - t62;
	t72 = t34 * t43;
	t70 = t37 * t45;
	t69 = t43 ^ 2 * t34;
	t48 = 0.1e1 / t50;
	t55 = 0.1e1 / t59;
	t68 = t48 * t55;
	t65 = t58 * t59;
	t44 = t57 * t66 + t63;
	t42 = 0.1e1 / t44 ^ 2;
	t61 = t58 ^ 2 * t73 * t42;
	t56 = 0.1e1 / t73;
	t49 = 0.1e1 / t50 ^ 2;
	t46 = t57 * t62 - t67;
	t41 = 0.1e1 / t44;
	t40 = 0.1e1 / (0.1e1 + t61);
	t38 = 0.1e1 / (t45 ^ 2 * t56 * t49 + 0.1e1);
	t33 = 0.1e1 / t35;
	t32 = (t45 * t48 * t56 * t57 + t60) * t38;
	t31 = 0.1e1 / (0.1e1 + t69);
	t30 = (-t45 * t49 * t51 + t46 * t48) * t55 * t38;
	t1 = [-t43 * t38 * t68, 0, t32, t30, 0, 0; (t45 * t33 - (-t36 + (-t68 * t70 + t36) * t38) * t69) * t31, 0, (-t32 * t70 * t72 + (t33 * t65 - (-t37 * t57 + (-t32 + t60) * t59 * t36) * t72) * t50) * t31, (t44 * t33 - (t37 * t59 * t51 + t36 * t46 + (-t36 * t64 + t70) * t30) * t72) * t31, 0, 0; (-t42 * t46 * t58 + t41 * t60) * t59 * t40, 0, (-t41 * t57 * t58 - t51 * t61) * t40, t43 * t42 * t40 * t65, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end