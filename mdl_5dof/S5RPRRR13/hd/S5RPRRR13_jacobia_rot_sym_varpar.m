% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRR13
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
%   Wie in S5RPRRR13_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPRRR13_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR13_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:02:48
	% EndTime: 2019-12-29 18:02:48
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:02:53
	% EndTime: 2019-12-29 18:02:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:02:43
	% EndTime: 2019-12-29 18:02:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:02:53
	% EndTime: 2019-12-29 18:02:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:02:48
	% EndTime: 2019-12-29 18:02:48
	% DurationCPUTime: 0.15s
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
	t1 = [t35 * t34 * t49, 0, t26, 0, 0; (-t22 * t45 + (-t32 * t35 * t37 * t48 + (-t34 + 0.1e1) * t42 * t31) * t56 * t55) * t21, 0, (t39 * t22 + (t43 * t53 + t32 * t42 + (-t32 * t45 - t53) * t26) * t55) * t40 * t21, 0, 0; ((t39 * t47 + t50) * t27 - (t39 * t46 - t51) * t54) * t25, 0, (t27 * t38 - t41 * t54) * t25 * t49, t44 * t25, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:02:48
	% EndTime: 2019-12-29 18:02:48
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (158->20), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->38)
	t54 = sin(qJ(1));
	t69 = t54 ^ 2;
	t53 = sin(qJ(3));
	t55 = cos(qJ(3));
	t56 = cos(qJ(1));
	t58 = t56 * t55;
	t45 = atan2(-t58, t53);
	t43 = sin(t45);
	t44 = cos(t45);
	t37 = -t43 * t58 + t44 * t53;
	t36 = 0.1e1 / t37 ^ 2;
	t68 = t36 * t55;
	t52 = qJ(4) + qJ(5);
	t47 = sin(t52);
	t60 = t56 * t47;
	t48 = cos(t52);
	t63 = t54 * t48;
	t42 = t53 * t63 + t60;
	t40 = 0.1e1 / t42 ^ 2;
	t59 = t56 * t48;
	t64 = t54 * t47;
	t41 = t53 * t64 - t59;
	t67 = t40 * t41;
	t66 = t43 * t53;
	t51 = t55 ^ 2;
	t65 = 0.1e1 / t53 ^ 2 * t51;
	t62 = t54 * t55;
	t46 = 0.1e1 / (t56 ^ 2 * t65 + 0.1e1);
	t61 = t56 * t46;
	t57 = t41 ^ 2 * t40 + 0.1e1;
	t49 = 0.1e1 / t53;
	t39 = 0.1e1 / t42;
	t38 = (0.1e1 + t65) * t61;
	t35 = 0.1e1 / t37;
	t34 = 0.1e1 / t57;
	t33 = 0.1e1 / (t69 * t51 * t36 + 0.1e1);
	t32 = t57 * t34;
	t1 = [t49 * t46 * t62, 0, t38, 0, 0; (-t35 * t58 + (-t44 * t49 * t51 * t61 + (-t46 + 0.1e1) * t55 * t43) * t69 * t68) * t33, 0, (t53 * t35 + (t56 * t66 + t44 * t55 + (-t44 * t58 - t66) * t38) * t68) * t54 * t33, 0, 0; ((t53 * t60 + t63) * t39 - (t53 * t59 - t64) * t67) * t34, 0, (t39 * t47 - t48 * t67) * t34 * t62, t32, t32;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end