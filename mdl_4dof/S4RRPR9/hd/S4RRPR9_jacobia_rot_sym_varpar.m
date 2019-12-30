% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S4RRPR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S4RRPR9_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR9_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_jacobia_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:46:42
	% EndTime: 2019-12-29 13:46:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:46:37
	% EndTime: 2019-12-29 13:46:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:46:44
	% EndTime: 2019-12-29 13:46:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:46:44
	% EndTime: 2019-12-29 13:46:44
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (99->20), mult. (197->54), div. (47->9), fcn. (297->9), ass. (0->32)
	t34 = cos(qJ(2));
	t32 = sin(qJ(2));
	t33 = sin(qJ(1));
	t39 = t33 * t32;
	t25 = atan2(-t39, -t34);
	t23 = sin(t25);
	t24 = cos(t25);
	t16 = -t23 * t39 - t24 * t34;
	t15 = 0.1e1 / t16 ^ 2;
	t35 = cos(qJ(1));
	t44 = t15 * t35 ^ 2;
	t30 = sin(pkin(7));
	t31 = cos(pkin(7));
	t36 = t35 * t31;
	t22 = t33 * t30 + t34 * t36;
	t20 = 0.1e1 / t22 ^ 2;
	t37 = t35 * t30;
	t21 = -t33 * t31 + t34 * t37;
	t43 = t20 * t21;
	t27 = t32 ^ 2;
	t42 = t27 / t34 ^ 2;
	t41 = t32 * t35;
	t26 = 0.1e1 / (t33 ^ 2 * t42 + 0.1e1);
	t40 = t33 * t26;
	t38 = t33 * t34;
	t28 = 0.1e1 / t34;
	t19 = 0.1e1 / t22;
	t18 = (0.1e1 + t42) * t40;
	t17 = 0.1e1 / (t21 ^ 2 * t20 + 0.1e1);
	t14 = 0.1e1 / t16;
	t13 = 0.1e1 / (t27 * t44 + 0.1e1);
	t1 = [t28 * t26 * t41, t18, 0, 0; (-t14 * t39 - (-t24 * t27 * t28 * t40 + (t26 - 0.1e1) * t32 * t23) * t32 * t44) * t13, (t34 * t14 - (-t23 * t38 + t24 * t32 + (t23 * t34 - t24 * t39) * t18) * t32 * t15) * t35 * t13, 0, 0; ((-t30 * t38 - t36) * t19 - (-t31 * t38 + t37) * t43) * t17, (-t19 * t30 + t31 * t43) * t17 * t41, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:46:44
	% EndTime: 2019-12-29 13:46:44
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (157->21), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->34)
	t46 = cos(qJ(2));
	t44 = sin(qJ(2));
	t45 = sin(qJ(1));
	t52 = t45 * t44;
	t36 = atan2(-t52, -t46);
	t34 = sin(t36);
	t35 = cos(t36);
	t28 = -t34 * t52 - t35 * t46;
	t27 = 0.1e1 / t28 ^ 2;
	t47 = cos(qJ(1));
	t57 = t27 * t47 ^ 2;
	t40 = pkin(7) + qJ(4);
	t38 = sin(t40);
	t39 = cos(t40);
	t49 = t47 * t39;
	t33 = t45 * t38 + t46 * t49;
	t31 = 0.1e1 / t33 ^ 2;
	t50 = t47 * t38;
	t32 = -t45 * t39 + t46 * t50;
	t56 = t31 * t32;
	t41 = t44 ^ 2;
	t55 = t41 / t46 ^ 2;
	t54 = t44 * t47;
	t37 = 0.1e1 / (t45 ^ 2 * t55 + 0.1e1);
	t53 = t45 * t37;
	t51 = t45 * t46;
	t48 = t32 ^ 2 * t31 + 0.1e1;
	t42 = 0.1e1 / t46;
	t30 = 0.1e1 / t33;
	t29 = (0.1e1 + t55) * t53;
	t26 = 0.1e1 / t28;
	t25 = 0.1e1 / t48;
	t24 = 0.1e1 / (t41 * t57 + 0.1e1);
	t1 = [t42 * t37 * t54, t29, 0, 0; (-t26 * t52 - (-t35 * t41 * t42 * t53 + (t37 - 0.1e1) * t44 * t34) * t44 * t57) * t24, (t46 * t26 - (-t34 * t51 + t35 * t44 + (t34 * t46 - t35 * t52) * t29) * t44 * t27) * t47 * t24, 0, 0; ((-t38 * t51 - t49) * t30 - (-t39 * t51 + t50) * t56) * t25, (-t30 * t38 + t39 * t56) * t25 * t54, 0, t48 * t25;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,4);
end