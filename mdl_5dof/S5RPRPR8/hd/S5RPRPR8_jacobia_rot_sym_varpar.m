% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRPR8
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
%   Wie in S5RPRPR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPRPR8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:51:51
	% EndTime: 2019-12-29 16:51:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:51:50
	% EndTime: 2019-12-29 16:51:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:51:51
	% EndTime: 2019-12-29 16:51:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:51:51
	% EndTime: 2019-12-29 16:51:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:51:56
	% EndTime: 2019-12-29 16:51:56
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (194->21), mult. (197->56), div. (47->9), fcn. (297->9), ass. (0->34)
	t30 = qJ(1) + pkin(8);
	t29 = cos(t30);
	t47 = t29 ^ 2;
	t37 = cos(qJ(3));
	t28 = sin(t30);
	t36 = sin(qJ(3));
	t42 = t28 * t36;
	t26 = atan2(-t42, -t37);
	t24 = sin(t26);
	t25 = cos(t26);
	t18 = -t24 * t42 - t25 * t37;
	t16 = 0.1e1 / t18 ^ 2;
	t46 = t16 * t36;
	t34 = sin(pkin(9));
	t35 = cos(pkin(9));
	t38 = t35 * t37;
	t23 = t28 * t34 + t29 * t38;
	t21 = 0.1e1 / t23 ^ 2;
	t39 = t34 * t37;
	t22 = -t28 * t35 + t29 * t39;
	t45 = t21 * t22;
	t44 = t24 * t37;
	t31 = t36 ^ 2;
	t40 = t31 / t37 ^ 2;
	t27 = 0.1e1 / (t28 ^ 2 * t40 + 0.1e1);
	t43 = t28 * t27;
	t41 = t29 * t36;
	t32 = 0.1e1 / t37;
	t20 = 0.1e1 / t23;
	t19 = (0.1e1 + t40) * t43;
	t17 = 0.1e1 / (t22 ^ 2 * t21 + 0.1e1);
	t15 = 0.1e1 / t18;
	t14 = 0.1e1 / (t47 * t31 * t16 + 0.1e1);
	t1 = [t32 * t27 * t41, 0, t19, 0, 0; (-t15 * t42 - (-t25 * t31 * t32 * t43 + (t27 - 0.1e1) * t36 * t24) * t47 * t46) * t14, 0, (t37 * t15 - (-t28 * t44 + t25 * t36 + (-t25 * t42 + t44) * t19) * t46) * t29 * t14, 0, 0; ((-t28 * t39 - t29 * t35) * t20 - (-t28 * t38 + t29 * t34) * t45) * t17, 0, (-t20 * t34 + t35 * t45) * t17 * t41, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:51:50
	% EndTime: 2019-12-29 16:51:51
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (266->22), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
	t44 = qJ(1) + pkin(8);
	t42 = cos(t44);
	t59 = t42 ^ 2;
	t49 = cos(qJ(3));
	t40 = sin(t44);
	t48 = sin(qJ(3));
	t55 = t40 * t48;
	t37 = atan2(-t55, -t49);
	t35 = sin(t37);
	t36 = cos(t37);
	t29 = -t35 * t55 - t36 * t49;
	t28 = 0.1e1 / t29 ^ 2;
	t58 = t28 * t48;
	t43 = pkin(9) + qJ(5);
	t39 = sin(t43);
	t41 = cos(t43);
	t52 = t42 * t49;
	t34 = t40 * t39 + t41 * t52;
	t32 = 0.1e1 / t34 ^ 2;
	t33 = t39 * t52 - t40 * t41;
	t57 = t32 * t33;
	t45 = t48 ^ 2;
	t51 = t45 / t49 ^ 2;
	t38 = 0.1e1 / (t40 ^ 2 * t51 + 0.1e1);
	t56 = t40 * t38;
	t54 = t40 * t49;
	t53 = t42 * t48;
	t50 = t33 ^ 2 * t32 + 0.1e1;
	t46 = 0.1e1 / t49;
	t31 = 0.1e1 / t34;
	t30 = (0.1e1 + t51) * t56;
	t27 = 0.1e1 / t29;
	t26 = 0.1e1 / t50;
	t25 = 0.1e1 / (t59 * t45 * t28 + 0.1e1);
	t1 = [t46 * t38 * t53, 0, t30, 0, 0; (-t27 * t55 - (-t36 * t45 * t46 * t56 + (t38 - 0.1e1) * t48 * t35) * t59 * t58) * t25, 0, (t49 * t27 - (-t35 * t54 + t36 * t48 + (t35 * t49 - t36 * t55) * t30) * t58) * t42 * t25, 0, 0; ((-t39 * t54 - t42 * t41) * t31 - (t42 * t39 - t41 * t54) * t57) * t26, 0, (-t31 * t39 + t41 * t57) * t26 * t53, 0, t50 * t26;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end