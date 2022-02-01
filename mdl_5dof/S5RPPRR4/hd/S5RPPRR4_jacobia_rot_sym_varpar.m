% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPRR4
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
%   Wie in S5RPPRR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPPRR4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
Ja_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (46->14), mult. (89->32), div. (20->9), fcn. (140->9), ass. (0->23)
	t28 = cos(pkin(8));
	t26 = sin(pkin(8));
	t29 = sin(qJ(1));
	t35 = t26 * t29;
	t21 = atan2(-t35, -t28);
	t19 = sin(t21);
	t20 = cos(t21);
	t14 = -t19 * t35 - t20 * t28;
	t30 = cos(qJ(1));
	t37 = 0.1e1 / t14 ^ 2 * t30 ^ 2;
	t23 = t26 ^ 2;
	t22 = 0.1e1 / (0.1e1 + t23 * t29 ^ 2 / t28 ^ 2);
	t36 = t22 / t28;
	t25 = sin(pkin(9));
	t34 = t29 * t25;
	t27 = cos(pkin(9));
	t33 = t29 * t27;
	t32 = t30 * t25;
	t31 = t30 * t27;
	t18 = t28 * t31 + t34;
	t17 = t28 * t32 - t33;
	t16 = 0.1e1 / t18 ^ 2;
	t1 = [t26 * t30 * t36, 0, 0, 0, 0; (-0.1e1 / t14 * t35 - (-t20 * t23 * t29 * t36 + (t22 - 0.1e1) * t26 * t19) * t26 * t37) / (t23 * t37 + 0.1e1), 0, 0, 0, 0; ((-t28 * t34 - t31) / t18 - (-t28 * t33 + t32) * t17 * t16) / (t17 ^ 2 * t16 + 0.1e1), 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (88->15), mult. (116->34), div. (25->9), fcn. (175->9), ass. (0->26)
	t34 = cos(pkin(8));
	t33 = sin(pkin(8));
	t35 = sin(qJ(1));
	t42 = t33 * t35;
	t26 = atan2(-t42, -t34);
	t24 = sin(t26);
	t25 = cos(t26);
	t20 = -t24 * t42 - t25 * t34;
	t36 = cos(qJ(1));
	t44 = 0.1e1 / t20 ^ 2 * t36 ^ 2;
	t30 = t33 ^ 2;
	t27 = 0.1e1 / (0.1e1 + t30 * t35 ^ 2 / t34 ^ 2);
	t43 = t27 / t34;
	t32 = pkin(9) + qJ(4);
	t28 = sin(t32);
	t41 = t35 * t28;
	t29 = cos(t32);
	t40 = t35 * t29;
	t39 = t36 * t28;
	t38 = t36 * t29;
	t23 = t34 * t38 + t41;
	t21 = 0.1e1 / t23 ^ 2;
	t22 = t34 * t39 - t40;
	t37 = t22 ^ 2 * t21 + 0.1e1;
	t18 = 0.1e1 / t37;
	t1 = [t33 * t36 * t43, 0, 0, 0, 0; (-0.1e1 / t20 * t42 - (-t25 * t30 * t35 * t43 + (t27 - 0.1e1) * t33 * t24) * t33 * t44) / (t30 * t44 + 0.1e1), 0, 0, 0, 0; ((-t34 * t41 - t38) / t23 - (-t34 * t40 + t39) * t22 * t21) * t18, 0, 0, t37 * t18, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.05s
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
	t47 = pkin(9) + qJ(4) + qJ(5);
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
	t1 = [t53 * t50 * t60, 0, 0, 0, 0; (-0.1e1 / t37 * t57 - (-t42 * t48 * t52 * t60 + (t44 - 0.1e1) * t50 * t41) * t50 * t61) / (t48 * t61 + 0.1e1), 0, 0, 0, 0; ((-t51 * t59 - t55) / t40 - (-t51 * t58 + t56) * t39 * t38) * t34, 0, 0, t33, t33;];
	Ja_rot = t1;
end