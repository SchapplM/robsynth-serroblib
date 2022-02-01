% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRPR5
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
%   Wie in S5RPRPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPRPR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
Ja_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (56->14), mult. (116->34), div. (25->9), fcn. (175->9), ass. (0->25)
	t31 = cos(pkin(8));
	t30 = sin(pkin(8));
	t33 = sin(qJ(1));
	t41 = t33 * t30;
	t26 = atan2(-t41, -t31);
	t24 = sin(t26);
	t25 = cos(t26);
	t19 = -t24 * t41 - t25 * t31;
	t35 = cos(qJ(1));
	t43 = 0.1e1 / t19 ^ 2 * t35 ^ 2;
	t28 = t30 ^ 2;
	t27 = 0.1e1 / (0.1e1 + t33 ^ 2 * t28 / t31 ^ 2);
	t42 = t27 / t31;
	t32 = sin(qJ(3));
	t40 = t33 * t32;
	t34 = cos(qJ(3));
	t39 = t33 * t34;
	t38 = t35 * t32;
	t37 = t35 * t34;
	t23 = t31 * t37 + t40;
	t21 = 0.1e1 / t23 ^ 2;
	t22 = t31 * t38 - t39;
	t36 = t22 ^ 2 * t21 + 0.1e1;
	t20 = 0.1e1 / t36;
	t1 = [t35 * t30 * t42, 0, 0, 0, 0; (-0.1e1 / t19 * t41 - (-t25 * t28 * t33 * t42 + (t27 - 0.1e1) * t30 * t24) * t30 * t43) / (t28 * t43 + 0.1e1), 0, 0, 0, 0; ((-t31 * t40 - t37) / t23 - (-t31 * t39 + t38) * t22 * t21) * t20, 0, t36 * t20, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (88->15), mult. (116->34), div. (25->9), fcn. (175->9), ass. (0->26)
	t37 = cos(pkin(8));
	t36 = sin(pkin(8));
	t38 = sin(qJ(1));
	t45 = t36 * t38;
	t29 = atan2(-t45, -t37);
	t27 = sin(t29);
	t28 = cos(t29);
	t23 = -t27 * t45 - t28 * t37;
	t39 = cos(qJ(1));
	t47 = 0.1e1 / t23 ^ 2 * t39 ^ 2;
	t33 = t36 ^ 2;
	t30 = 0.1e1 / (0.1e1 + t33 * t38 ^ 2 / t37 ^ 2);
	t46 = t30 / t37;
	t35 = qJ(3) + pkin(9);
	t31 = sin(t35);
	t44 = t38 * t31;
	t32 = cos(t35);
	t43 = t38 * t32;
	t42 = t39 * t31;
	t41 = t39 * t32;
	t26 = t37 * t41 + t44;
	t24 = 0.1e1 / t26 ^ 2;
	t25 = t37 * t42 - t43;
	t40 = t25 ^ 2 * t24 + 0.1e1;
	t21 = 0.1e1 / t40;
	t1 = [t36 * t39 * t46, 0, 0, 0, 0; (-0.1e1 / t23 * t45 - (-t28 * t33 * t38 * t46 + (t30 - 0.1e1) * t36 * t27) * t36 * t47) / (t33 * t47 + 0.1e1), 0, 0, 0, 0; ((-t37 * t44 - t41) / t26 - (-t37 * t43 + t42) * t25 * t24) * t21, 0, t40 * t21, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (158->15), mult. (143->34), div. (30->9), fcn. (210->9), ass. (0->27)
	t51 = cos(pkin(8));
	t50 = sin(pkin(8));
	t52 = sin(qJ(1));
	t57 = t52 * t50;
	t43 = atan2(-t57, -t51);
	t41 = sin(t43);
	t42 = cos(t43);
	t37 = -t41 * t57 - t42 * t51;
	t53 = cos(qJ(1));
	t61 = 0.1e1 / t37 ^ 2 * t53 ^ 2;
	t48 = t50 ^ 2;
	t44 = 0.1e1 / (0.1e1 + t52 ^ 2 * t48 / t51 ^ 2);
	t60 = t44 / t51;
	t47 = qJ(3) + pkin(9) + qJ(5);
	t45 = sin(t47);
	t59 = t52 * t45;
	t46 = cos(t47);
	t58 = t52 * t46;
	t56 = t53 * t45;
	t55 = t53 * t46;
	t40 = t51 * t55 + t59;
	t38 = 0.1e1 / t40 ^ 2;
	t39 = t51 * t56 - t58;
	t54 = t39 ^ 2 * t38 + 0.1e1;
	t34 = 0.1e1 / t54;
	t33 = t54 * t34;
	t1 = [t53 * t50 * t60, 0, 0, 0, 0; (-0.1e1 / t37 * t57 - (-t42 * t48 * t52 * t60 + (t44 - 0.1e1) * t50 * t41) * t50 * t61) / (t48 * t61 + 0.1e1), 0, 0, 0, 0; ((-t51 * t59 - t55) / t40 - (-t51 * t58 + t56) * t39 * t38) * t34, 0, t33, 0, t33;];
	Ja_rot = t1;
end