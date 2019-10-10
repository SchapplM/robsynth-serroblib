% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR2
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
%   Wie in S6RPRRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (218->21), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
	t36 = qJ(1) + pkin(10);
	t35 = cos(t36);
	t54 = t35 ^ 2;
	t43 = cos(qJ(3));
	t34 = sin(t36);
	t41 = sin(qJ(3));
	t49 = t34 * t41;
	t32 = atan2(-t49, -t43);
	t30 = sin(t32);
	t31 = cos(t32);
	t23 = -t30 * t49 - t31 * t43;
	t22 = 0.1e1 / t23 ^ 2;
	t53 = t22 * t41;
	t40 = sin(qJ(4));
	t42 = cos(qJ(4));
	t45 = t42 * t43;
	t29 = t34 * t40 + t35 * t45;
	t27 = 0.1e1 / t29 ^ 2;
	t46 = t40 * t43;
	t28 = -t34 * t42 + t35 * t46;
	t52 = t27 * t28;
	t51 = t30 * t43;
	t37 = t41 ^ 2;
	t47 = t37 / t43 ^ 2;
	t33 = 0.1e1 / (t34 ^ 2 * t47 + 0.1e1);
	t50 = t34 * t33;
	t48 = t35 * t41;
	t44 = t28 ^ 2 * t27 + 0.1e1;
	t38 = 0.1e1 / t43;
	t26 = 0.1e1 / t29;
	t25 = (0.1e1 + t47) * t50;
	t24 = 0.1e1 / t44;
	t21 = 0.1e1 / t23;
	t20 = 0.1e1 / (t54 * t37 * t22 + 0.1e1);
	t1 = [t38 * t33 * t48, 0, t25, 0, 0, 0; (-t21 * t49 - (-t31 * t37 * t38 * t50 + (t33 - 0.1e1) * t41 * t30) * t54 * t53) * t20, 0, (t43 * t21 - (-t34 * t51 + t31 * t41 + (-t31 * t49 + t51) * t25) * t53) * t35 * t20, 0, 0, 0; ((-t34 * t46 - t35 * t42) * t26 - (-t34 * t45 + t35 * t40) * t52) * t24, 0, (-t26 * t40 + t42 * t52) * t24 * t48, t44 * t24, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (266->22), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
	t47 = qJ(1) + pkin(10);
	t45 = cos(t47);
	t62 = t45 ^ 2;
	t52 = cos(qJ(3));
	t43 = sin(t47);
	t51 = sin(qJ(3));
	t58 = t43 * t51;
	t40 = atan2(-t58, -t52);
	t38 = sin(t40);
	t39 = cos(t40);
	t32 = -t38 * t58 - t39 * t52;
	t31 = 0.1e1 / t32 ^ 2;
	t61 = t31 * t51;
	t46 = qJ(4) + pkin(11);
	t42 = sin(t46);
	t44 = cos(t46);
	t55 = t45 * t52;
	t37 = t43 * t42 + t44 * t55;
	t35 = 0.1e1 / t37 ^ 2;
	t36 = t42 * t55 - t43 * t44;
	t60 = t35 * t36;
	t48 = t51 ^ 2;
	t54 = t48 / t52 ^ 2;
	t41 = 0.1e1 / (t43 ^ 2 * t54 + 0.1e1);
	t59 = t43 * t41;
	t57 = t43 * t52;
	t56 = t45 * t51;
	t53 = t36 ^ 2 * t35 + 0.1e1;
	t49 = 0.1e1 / t52;
	t34 = 0.1e1 / t37;
	t33 = (0.1e1 + t54) * t59;
	t30 = 0.1e1 / t32;
	t29 = 0.1e1 / t53;
	t28 = 0.1e1 / (t62 * t48 * t31 + 0.1e1);
	t1 = [t49 * t41 * t56, 0, t33, 0, 0, 0; (-t30 * t58 - (-t39 * t48 * t49 * t59 + (t41 - 0.1e1) * t51 * t38) * t62 * t61) * t28, 0, (t52 * t30 - (-t38 * t57 + t39 * t51 + (t38 * t52 - t39 * t58) * t33) * t61) * t45 * t28, 0, 0, 0; ((-t42 * t57 - t45 * t44) * t34 - (t45 * t42 - t44 * t57) * t60) * t29, 0, (-t34 * t42 + t44 * t60) * t29 * t56, t53 * t29, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:25:14
	% EndTime: 2019-10-10 01:25:14
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (366->22), mult. (251->57), div. (57->9), fcn. (367->9), ass. (0->36)
	t56 = qJ(1) + pkin(10);
	t54 = cos(t56);
	t71 = t54 ^ 2;
	t61 = cos(qJ(3));
	t53 = sin(t56);
	t60 = sin(qJ(3));
	t67 = t53 * t60;
	t49 = atan2(-t67, -t61);
	t47 = sin(t49);
	t48 = cos(t49);
	t41 = -t47 * t67 - t48 * t61;
	t40 = 0.1e1 / t41 ^ 2;
	t70 = t40 * t60;
	t55 = qJ(4) + pkin(11) + qJ(6);
	t51 = sin(t55);
	t52 = cos(t55);
	t64 = t54 * t61;
	t46 = t53 * t51 + t52 * t64;
	t44 = 0.1e1 / t46 ^ 2;
	t45 = t51 * t64 - t53 * t52;
	t69 = t44 * t45;
	t57 = t60 ^ 2;
	t63 = t57 / t61 ^ 2;
	t50 = 0.1e1 / (t53 ^ 2 * t63 + 0.1e1);
	t68 = t53 * t50;
	t66 = t53 * t61;
	t65 = t54 * t60;
	t62 = t45 ^ 2 * t44 + 0.1e1;
	t58 = 0.1e1 / t61;
	t43 = 0.1e1 / t46;
	t42 = (0.1e1 + t63) * t68;
	t39 = 0.1e1 / t41;
	t38 = 0.1e1 / (t71 * t57 * t40 + 0.1e1);
	t37 = 0.1e1 / t62;
	t36 = t62 * t37;
	t1 = [t58 * t50 * t65, 0, t42, 0, 0, 0; (-t39 * t67 - (-t48 * t57 * t58 * t68 + (t50 - 0.1e1) * t60 * t47) * t71 * t70) * t38, 0, (t61 * t39 - (-t47 * t66 + t48 * t60 + (t47 * t61 - t48 * t67) * t42) * t70) * t54 * t38, 0, 0, 0; ((-t51 * t66 - t54 * t52) * t43 - (t54 * t51 - t52 * t66) * t69) * t37, 0, (-t43 * t51 + t52 * t69) * t37 * t65, t36, 0, t36;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end