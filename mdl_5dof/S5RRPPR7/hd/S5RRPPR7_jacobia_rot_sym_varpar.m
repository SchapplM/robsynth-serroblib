% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPPR7
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
%   Wie in S5RRPPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRPPR7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (192->17), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->28)
	t36 = cos(qJ(1));
	t34 = t36 ^ 2;
	t32 = qJ(2) + pkin(8);
	t31 = cos(t32);
	t30 = sin(t32);
	t35 = sin(qJ(1));
	t39 = t35 * t30;
	t24 = atan2(-t39, -t31);
	t22 = sin(t24);
	t23 = cos(t24);
	t20 = -t22 * t39 - t23 * t31;
	t19 = 0.1e1 / t20 ^ 2;
	t45 = t19 * t30;
	t44 = t22 * t31;
	t27 = t30 ^ 2;
	t37 = t31 ^ 2;
	t43 = t27 / t37;
	t42 = t30 * t36;
	t38 = t35 ^ 2;
	t41 = 0.1e1 / t38 * t34;
	t25 = 0.1e1 / (t38 * t43 + 0.1e1);
	t40 = t35 * t25;
	t28 = 0.1e1 / t31;
	t26 = 0.1e1 / (t37 * t41 + 0.1e1);
	t21 = (0.1e1 + t43) * t40;
	t18 = 0.1e1 / t20;
	t17 = 0.1e1 / (t34 * t27 * t19 + 0.1e1);
	t1 = [t28 * t25 * t42, t21, 0, 0, 0; (-t18 * t39 - (-t23 * t27 * t28 * t40 + (t25 - 0.1e1) * t30 * t22) * t34 * t45) * t17, (t31 * t18 - (-t35 * t44 + t23 * t30 + (-t23 * t39 + t44) * t21) * t45) * t36 * t17, 0, 0, 0; (-0.1e1 - t41) * t31 * t26, -0.1e1 / t35 * t26 * t42, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (216->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
	t47 = qJ(2) + pkin(8);
	t45 = sin(t47);
	t46 = cos(t47);
	t49 = sin(qJ(1));
	t57 = t49 * t46;
	t40 = atan2(-t57, t45);
	t38 = sin(t40);
	t39 = cos(t40);
	t31 = -t38 * t57 + t39 * t45;
	t30 = 0.1e1 / t31 ^ 2;
	t51 = cos(qJ(1));
	t63 = t30 * t51 ^ 2;
	t48 = sin(qJ(5));
	t54 = t51 * t48;
	t50 = cos(qJ(5));
	t55 = t49 * t50;
	t37 = t45 * t54 + t55;
	t35 = 0.1e1 / t37 ^ 2;
	t53 = t51 * t50;
	t56 = t49 * t48;
	t36 = -t45 * t53 + t56;
	t62 = t35 * t36;
	t61 = t38 * t45;
	t44 = t46 ^ 2;
	t60 = 0.1e1 / t45 ^ 2 * t44;
	t59 = t46 * t51;
	t41 = 0.1e1 / (t49 ^ 2 * t60 + 0.1e1);
	t58 = t49 * t41;
	t52 = t36 ^ 2 * t35 + 0.1e1;
	t42 = 0.1e1 / t45;
	t34 = 0.1e1 / t37;
	t33 = 0.1e1 / t52;
	t32 = (0.1e1 + t60) * t58;
	t29 = 0.1e1 / t31;
	t28 = 0.1e1 / (t44 * t63 + 0.1e1);
	t1 = [-t42 * t41 * t59, t32, 0, 0, 0; (-t29 * t57 - (t39 * t42 * t44 * t58 + (t41 - 0.1e1) * t46 * t38) * t46 * t63) * t28, (-t45 * t29 - (t49 * t61 + t39 * t46 + (-t39 * t57 - t61) * t32) * t46 * t30) * t51 * t28, 0, 0, 0; ((t45 * t55 + t54) * t34 - (-t45 * t56 + t53) * t62) * t33, (-t34 * t50 - t48 * t62) * t33 * t59, 0, 0, t52 * t33;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end