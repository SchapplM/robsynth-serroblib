% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPRR3
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
%   Wie in S5PRPRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRPRR3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:48:12
	% EndTime: 2019-12-05 15:48:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:48:12
	% EndTime: 2019-12-05 15:48:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:48:12
	% EndTime: 2019-12-05 15:48:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:48:12
	% EndTime: 2019-12-05 15:48:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:48:12
	% EndTime: 2019-12-05 15:48:12
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (137->15), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->24)
	t32 = qJ(2) + pkin(9);
	t31 = cos(t32);
	t30 = sin(t32);
	t33 = sin(pkin(8));
	t38 = t33 * t30;
	t26 = atan2(-t38, -t31);
	t24 = sin(t26);
	t41 = t24 * t31;
	t28 = t30 ^ 2;
	t40 = t28 / t31 ^ 2;
	t34 = cos(pkin(8));
	t39 = t31 * t34;
	t35 = sin(qJ(4));
	t36 = cos(qJ(4));
	t23 = t33 * t35 + t36 * t39;
	t21 = 0.1e1 / t23 ^ 2;
	t22 = -t33 * t36 + t35 * t39;
	t37 = t22 ^ 2 * t21 + 0.1e1;
	t25 = cos(t26);
	t20 = 0.1e1 / t37;
	t19 = (0.1e1 + t40) * t33 / (t33 ^ 2 * t40 + 0.1e1);
	t18 = -t24 * t38 - t25 * t31;
	t17 = 0.1e1 / t18 ^ 2;
	t1 = [0, t19, 0, 0, 0; 0, (t31 / t18 - (-t33 * t41 + t25 * t30 + (-t25 * t38 + t41) * t19) * t30 * t17) * t34 / (t34 ^ 2 * t28 * t17 + 0.1e1), 0, 0, 0; 0, (-t35 / t23 + t36 * t22 * t21) * t34 * t30 * t20, 0, t37 * t20, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:48:12
	% EndTime: 2019-12-05 15:48:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (198->16), mult. (162->34), div. (37->8), fcn. (227->9), ass. (0->26)
	t48 = qJ(2) + pkin(9);
	t45 = cos(t48);
	t44 = sin(t48);
	t50 = sin(pkin(8));
	t53 = t50 * t44;
	t40 = atan2(-t53, -t45);
	t38 = sin(t40);
	t56 = t38 * t45;
	t42 = t44 ^ 2;
	t55 = t42 / t45 ^ 2;
	t51 = cos(pkin(8));
	t54 = t45 * t51;
	t49 = qJ(4) + qJ(5);
	t46 = sin(t49);
	t47 = cos(t49);
	t37 = t50 * t46 + t47 * t54;
	t35 = 0.1e1 / t37 ^ 2;
	t36 = t46 * t54 - t50 * t47;
	t52 = t36 ^ 2 * t35 + 0.1e1;
	t39 = cos(t40);
	t34 = 0.1e1 / t52;
	t33 = (0.1e1 + t55) * t50 / (t50 ^ 2 * t55 + 0.1e1);
	t32 = -t38 * t53 - t39 * t45;
	t31 = 0.1e1 / t32 ^ 2;
	t29 = t52 * t34;
	t1 = [0, t33, 0, 0, 0; 0, (t45 / t32 - (-t50 * t56 + t39 * t44 + (-t39 * t53 + t56) * t33) * t44 * t31) * t51 / (t51 ^ 2 * t42 * t31 + 0.1e1), 0, 0, 0; 0, (-t46 / t37 + t47 * t36 * t35) * t51 * t44 * t34, 0, t29, t29;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end