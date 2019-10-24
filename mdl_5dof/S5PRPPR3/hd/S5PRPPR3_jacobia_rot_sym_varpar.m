% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPPR3
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
%   Wie in S5PRPPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:22
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRPPR3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:52
	% EndTime: 2019-10-24 10:22:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:52
	% EndTime: 2019-10-24 10:22:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:52
	% EndTime: 2019-10-24 10:22:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:52
	% EndTime: 2019-10-24 10:22:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:52
	% EndTime: 2019-10-24 10:22:52
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (105->13), mult. (82->22), div. (24->8), fcn. (122->7), ass. (0->19)
	t23 = qJ(2) + pkin(8);
	t21 = cos(t23);
	t20 = sin(t23);
	t24 = sin(pkin(7));
	t28 = t24 * t20;
	t16 = atan2(-t28, -t21);
	t14 = sin(t16);
	t30 = t14 * t21;
	t18 = t20 ^ 2;
	t26 = t21 ^ 2;
	t29 = t18 / t26;
	t27 = t24 ^ 2;
	t25 = cos(pkin(7));
	t22 = t25 ^ 2;
	t15 = cos(t16);
	t13 = (0.1e1 + t29) * t24 / (t27 * t29 + 0.1e1);
	t12 = -t14 * t28 - t15 * t21;
	t11 = 0.1e1 / t12 ^ 2;
	t1 = [0, t13, 0, 0, 0; 0, (t21 / t12 - (-t24 * t30 + t15 * t20 + (-t15 * t28 + t30) * t13) * t20 * t11) * t25 / (t22 * t18 * t11 + 0.1e1), 0, 0, 0; 0, -t25 * t20 / t24 / (0.1e1 + t22 * t26 / t27), 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:52
	% EndTime: 2019-10-24 10:22:52
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (125->14), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->24)
	t32 = qJ(2) + pkin(8);
	t30 = sin(t32);
	t31 = cos(t32);
	t33 = sin(pkin(7));
	t38 = t33 * t31;
	t26 = atan2(-t38, t30);
	t24 = sin(t26);
	t41 = t24 * t30;
	t29 = t31 ^ 2;
	t40 = 0.1e1 / t30 ^ 2 * t29;
	t34 = cos(pkin(7));
	t39 = t30 * t34;
	t35 = sin(qJ(5));
	t36 = cos(qJ(5));
	t23 = t33 * t36 + t35 * t39;
	t21 = 0.1e1 / t23 ^ 2;
	t22 = t33 * t35 - t36 * t39;
	t37 = t22 ^ 2 * t21 + 0.1e1;
	t25 = cos(t26);
	t20 = 0.1e1 / t37;
	t19 = (0.1e1 + t40) * t33 / (t33 ^ 2 * t40 + 0.1e1);
	t18 = -t24 * t38 + t25 * t30;
	t17 = 0.1e1 / t18 ^ 2;
	t1 = [0, t19, 0, 0, 0; 0, (-t30 / t18 - (t33 * t41 + t25 * t31 + (-t25 * t38 - t41) * t19) * t31 * t17) * t34 / (t34 ^ 2 * t29 * t17 + 0.1e1), 0, 0, 0; 0, (-t36 / t23 - t35 * t22 * t21) * t34 * t31 * t20, 0, 0, t37 * t20;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end