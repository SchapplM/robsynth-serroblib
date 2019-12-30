% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPPR5
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
%   Wie in S5PRPPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRPPR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:24
	% EndTime: 2019-12-29 15:31:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:30
	% EndTime: 2019-12-29 15:31:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:19
	% EndTime: 2019-12-29 15:31:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:24
	% EndTime: 2019-12-29 15:31:24
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (45->12), mult. (82->23), div. (26->8), fcn. (123->7), ass. (0->18)
	t22 = cos(qJ(2));
	t19 = sin(pkin(7));
	t21 = sin(qJ(2));
	t24 = t19 * t21;
	t14 = atan2(-t24, -t22);
	t12 = sin(t14);
	t26 = t12 * t22;
	t17 = t21 ^ 2;
	t18 = 0.1e1 / t22 ^ 2;
	t25 = t17 * t18;
	t20 = cos(pkin(7));
	t23 = t20 ^ 2;
	t16 = t19 ^ 2;
	t13 = cos(t14);
	t11 = (0.1e1 + t25) * t19 / (t16 * t25 + 0.1e1);
	t10 = -t12 * t24 - t13 * t22;
	t9 = 0.1e1 / t10 ^ 2;
	t1 = [0, t11, 0, 0, 0; 0, (t22 / t10 - (-t19 * t26 + t13 * t21 + (-t13 * t24 + t26) * t11) * t21 * t9) / (t23 * t17 * t9 + 0.1e1) * t20, 0, 0, 0; 0, -0.1e1 / t20 * t18 / (0.1e1 + t16 / t23 * t18) * t24, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:25
	% EndTime: 2019-12-29 15:31:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:31:24
	% EndTime: 2019-12-29 15:31:25
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (122->15), mult. (346->37), div. (35->9), fcn. (500->11), ass. (0->27)
	t52 = sin(pkin(8));
	t54 = cos(pkin(8));
	t57 = sin(qJ(2));
	t59 = cos(qJ(2));
	t47 = t57 * t52 + t59 * t54;
	t48 = t59 * t52 - t57 * t54;
	t55 = cos(pkin(7));
	t44 = t47 * t55;
	t53 = sin(pkin(7));
	t56 = sin(qJ(5));
	t58 = cos(qJ(5));
	t40 = t44 * t58 - t53 * t56;
	t38 = 0.1e1 / t40 ^ 2;
	t39 = t44 * t56 + t53 * t58;
	t60 = t39 ^ 2 * t38 + 0.1e1;
	t46 = 0.1e1 / t47 ^ 2;
	t45 = t48 * t55;
	t42 = t48 * t53;
	t41 = t47 * t53;
	t37 = atan2(-t42, t47);
	t35 = cos(t37);
	t34 = sin(t37);
	t33 = 0.1e1 / t60;
	t32 = -t34 * t42 + t35 * t47;
	t31 = 0.1e1 / t32 ^ 2;
	t29 = (t41 / t47 + t48 * t42 * t46) / (t42 ^ 2 * t46 + 0.1e1);
	t1 = [0, t29, 0, 0, 0; 0, (-t44 / t32 - (t34 * t41 + t35 * t48 + (-t34 * t47 - t35 * t42) * t29) * t45 * t31) / (t45 ^ 2 * t31 + 0.1e1), 0, 0, 0; 0, (t56 / t40 - t58 * t39 * t38) * t45 * t33, 0, 0, t60 * t33;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end