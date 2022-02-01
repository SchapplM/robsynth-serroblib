% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRP4
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
%   Wie in S5RPRRP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPRRP4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
Ja_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
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
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (112->15), mult. (143->34), div. (30->9), fcn. (210->9), ass. (0->27)
	t46 = cos(pkin(8));
	t45 = sin(pkin(8));
	t47 = sin(qJ(1));
	t52 = t47 * t45;
	t38 = atan2(-t52, -t46);
	t36 = sin(t38);
	t37 = cos(t38);
	t32 = -t36 * t52 - t37 * t46;
	t48 = cos(qJ(1));
	t56 = 0.1e1 / t32 ^ 2 * t48 ^ 2;
	t42 = t45 ^ 2;
	t39 = 0.1e1 / (0.1e1 + t47 ^ 2 * t42 / t46 ^ 2);
	t55 = t39 / t46;
	t44 = qJ(3) + qJ(4);
	t40 = sin(t44);
	t54 = t47 * t40;
	t41 = cos(t44);
	t53 = t47 * t41;
	t51 = t48 * t40;
	t50 = t48 * t41;
	t35 = t46 * t50 + t54;
	t33 = 0.1e1 / t35 ^ 2;
	t34 = t46 * t51 - t53;
	t49 = t34 ^ 2 * t33 + 0.1e1;
	t30 = 0.1e1 / t49;
	t28 = t49 * t30;
	t1 = [t48 * t45 * t55, 0, 0, 0, 0; (-0.1e1 / t32 * t52 - (-t37 * t42 * t47 * t55 + (t39 - 0.1e1) * t45 * t36) * t45 * t56) / (t42 * t56 + 0.1e1), 0, 0, 0, 0; ((-t46 * t54 - t50) / t35 - (-t46 * t53 + t51) * t34 * t33) * t30, 0, t28, t28, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (112->15), mult. (143->34), div. (30->9), fcn. (210->9), ass. (0->27)
	t50 = cos(pkin(8));
	t49 = sin(pkin(8));
	t51 = sin(qJ(1));
	t56 = t51 * t49;
	t42 = atan2(-t56, -t50);
	t40 = sin(t42);
	t41 = cos(t42);
	t36 = -t40 * t56 - t41 * t50;
	t52 = cos(qJ(1));
	t60 = 0.1e1 / t36 ^ 2 * t52 ^ 2;
	t46 = t49 ^ 2;
	t43 = 0.1e1 / (0.1e1 + t51 ^ 2 * t46 / t50 ^ 2);
	t59 = t43 / t50;
	t48 = qJ(3) + qJ(4);
	t44 = sin(t48);
	t58 = t51 * t44;
	t45 = cos(t48);
	t57 = t51 * t45;
	t55 = t52 * t44;
	t54 = t52 * t45;
	t39 = t50 * t54 + t58;
	t37 = 0.1e1 / t39 ^ 2;
	t38 = t50 * t55 - t57;
	t53 = t38 ^ 2 * t37 + 0.1e1;
	t34 = 0.1e1 / t53;
	t32 = t53 * t34;
	t1 = [t52 * t49 * t59, 0, 0, 0, 0; (-0.1e1 / t36 * t56 - (-t41 * t46 * t51 * t59 + (t43 - 0.1e1) * t49 * t40) * t49 * t60) / (t46 * t60 + 0.1e1), 0, 0, 0, 0; ((-t50 * t58 - t54) / t39 - (-t50 * t57 + t55) * t38 * t37) * t34, 0, t32, t32, 0;];
	Ja_rot = t1;
end