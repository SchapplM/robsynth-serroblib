% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR8
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
%   Wie in S6RPRRPR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:38
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:38
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
	t41 = qJ(4) + pkin(10);
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
	% StartTime: 2019-10-10 01:35:37
	% EndTime: 2019-10-10 01:35:38
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (220->20), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->38)
	t59 = sin(qJ(1));
	t74 = t59 ^ 2;
	t58 = sin(qJ(3));
	t60 = cos(qJ(3));
	t61 = cos(qJ(1));
	t63 = t61 * t60;
	t50 = atan2(-t63, t58);
	t48 = sin(t50);
	t49 = cos(t50);
	t42 = -t48 * t63 + t49 * t58;
	t41 = 0.1e1 / t42 ^ 2;
	t73 = t41 * t60;
	t54 = qJ(4) + pkin(10) + qJ(6);
	t52 = sin(t54);
	t65 = t61 * t52;
	t53 = cos(t54);
	t68 = t59 * t53;
	t47 = t58 * t68 + t65;
	t45 = 0.1e1 / t47 ^ 2;
	t64 = t61 * t53;
	t69 = t59 * t52;
	t46 = t58 * t69 - t64;
	t72 = t45 * t46;
	t71 = t48 * t58;
	t57 = t60 ^ 2;
	t70 = 0.1e1 / t58 ^ 2 * t57;
	t67 = t59 * t60;
	t51 = 0.1e1 / (t61 ^ 2 * t70 + 0.1e1);
	t66 = t61 * t51;
	t62 = t46 ^ 2 * t45 + 0.1e1;
	t55 = 0.1e1 / t58;
	t44 = 0.1e1 / t47;
	t43 = (0.1e1 + t70) * t66;
	t40 = 0.1e1 / t42;
	t39 = 0.1e1 / (t74 * t57 * t41 + 0.1e1);
	t38 = 0.1e1 / t62;
	t37 = t62 * t38;
	t1 = [t55 * t51 * t67, 0, t43, 0, 0, 0; (-t40 * t63 + (-t49 * t55 * t57 * t66 + (-t51 + 0.1e1) * t60 * t48) * t74 * t73) * t39, 0, (t58 * t40 + (t61 * t71 + t49 * t60 + (-t49 * t63 - t71) * t43) * t73) * t59 * t39, 0, 0, 0; ((t58 * t65 + t68) * t44 - (t58 * t64 - t69) * t72) * t38, 0, (t44 * t52 - t53 * t72) * t38 * t67, t37, 0, t37;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end