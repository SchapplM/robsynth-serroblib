% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4PRRP6
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
%   Wie in S4PRRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S4PRRP6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRP6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_jacobia_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:22:13
	% EndTime: 2019-12-29 12:22:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:22:07
	% EndTime: 2019-12-29 12:22:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:22:07
	% EndTime: 2019-12-29 12:22:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:22:13
	% EndTime: 2019-12-29 12:22:13
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (63->14), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->23)
	t31 = cos(qJ(2));
	t26 = sin(pkin(6));
	t29 = sin(qJ(2));
	t34 = t26 * t29;
	t22 = atan2(-t34, -t31);
	t20 = sin(t22);
	t36 = t20 * t31;
	t24 = t29 ^ 2;
	t35 = t24 / t31 ^ 2;
	t27 = cos(pkin(6));
	t33 = t27 * t31;
	t28 = sin(qJ(3));
	t30 = cos(qJ(3));
	t19 = t26 * t28 + t30 * t33;
	t17 = 0.1e1 / t19 ^ 2;
	t18 = -t26 * t30 + t28 * t33;
	t32 = t18 ^ 2 * t17 + 0.1e1;
	t21 = cos(t22);
	t16 = (0.1e1 + t35) * t26 / (t26 ^ 2 * t35 + 0.1e1);
	t15 = 0.1e1 / t32;
	t14 = -t20 * t34 - t21 * t31;
	t13 = 0.1e1 / t14 ^ 2;
	t1 = [0, t16, 0, 0; 0, (t31 / t14 - (-t26 * t36 + t21 * t29 + (-t21 * t34 + t36) * t16) * t29 * t13) * t27 / (t27 ^ 2 * t24 * t13 + 0.1e1), 0, 0; 0, (-t28 / t19 + t30 * t18 * t17) * t29 * t27 * t15, t32 * t15, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:22:13
	% EndTime: 2019-12-29 12:22:13
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (112->21), mult. (359->57), div. (75->11), fcn. (540->9), ass. (0->34)
	t45 = sin(qJ(2));
	t56 = t45 ^ 2;
	t42 = sin(pkin(6));
	t43 = cos(pkin(6));
	t46 = cos(qJ(3));
	t44 = sin(qJ(3));
	t47 = cos(qJ(2));
	t51 = t44 * t47;
	t31 = t42 * t51 + t43 * t46;
	t50 = t45 * t44;
	t29 = atan2(-t31, t50);
	t25 = sin(t29);
	t26 = cos(t29);
	t24 = -t25 * t31 + t26 * t50;
	t23 = 0.1e1 / t24 ^ 2;
	t34 = -t42 * t46 + t43 * t51;
	t55 = t23 * t34;
	t53 = t26 * t31;
	t52 = t43 * t45;
	t49 = t46 * t47;
	t35 = t42 * t44 + t43 * t49;
	t30 = 0.1e1 / t35 ^ 2;
	t48 = t43 ^ 2 * t56 * t30;
	t41 = 0.1e1 / t56;
	t38 = 0.1e1 / t44 ^ 2;
	t37 = 0.1e1 / t44;
	t33 = t42 * t49 - t43 * t44;
	t28 = 0.1e1 / (t31 ^ 2 * t41 * t38 + 0.1e1);
	t27 = 0.1e1 / (0.1e1 + t48);
	t22 = 0.1e1 / t24;
	t21 = (t31 * t37 * t41 * t47 + t42) * t28;
	t20 = 0.1e1 / (t34 ^ 2 * t23 + 0.1e1);
	t19 = (t31 * t38 * t46 - t33 * t37) / t45 * t28;
	t1 = [0, t21, t19, 0; 0, (t21 * t53 * t55 + (-t22 * t52 - (t26 * t47 + (-t21 + t42) * t45 * t25) * t55) * t44) * t20, (t35 * t22 - (t26 * t45 * t46 - t25 * t33 + (-t25 * t50 - t53) * t19) * t55) * t20, 0; 0, (-t43 * t47 / t35 - t46 * t48) * t27, -t34 * t30 * t27 * t52, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,4);
end