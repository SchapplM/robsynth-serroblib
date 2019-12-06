% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRPR2
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
%   Wie in S5PPRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PPRPR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (105->13), mult. (82->22), div. (24->8), fcn. (122->7), ass. (0->19)
	t22 = pkin(8) + qJ(3);
	t20 = cos(t22);
	t19 = sin(t22);
	t23 = sin(pkin(7));
	t27 = t23 * t19;
	t15 = atan2(-t27, -t20);
	t13 = sin(t15);
	t29 = t13 * t20;
	t17 = t19 ^ 2;
	t25 = t20 ^ 2;
	t28 = t17 / t25;
	t26 = t23 ^ 2;
	t24 = cos(pkin(7));
	t21 = t24 ^ 2;
	t14 = cos(t15);
	t12 = (0.1e1 + t28) * t23 / (t26 * t28 + 0.1e1);
	t11 = -t13 * t27 - t14 * t20;
	t10 = 0.1e1 / t11 ^ 2;
	t1 = [0, 0, t12, 0, 0; 0, 0, (t20 / t11 - (-t23 * t29 + t14 * t19 + (-t14 * t27 + t29) * t12) * t19 * t10) / (t21 * t17 * t10 + 0.1e1) * t24, 0, 0; 0, 0, -t24 * t19 / t23 / (0.1e1 + t21 * t25 / t26), 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:03:38
	% EndTime: 2019-12-05 15:03:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (125->14), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->24)
	t31 = pkin(8) + qJ(3);
	t29 = sin(t31);
	t30 = cos(t31);
	t32 = sin(pkin(7));
	t37 = t32 * t30;
	t25 = atan2(-t37, t29);
	t23 = sin(t25);
	t40 = t23 * t29;
	t28 = t30 ^ 2;
	t39 = 0.1e1 / t29 ^ 2 * t28;
	t33 = cos(pkin(7));
	t38 = t29 * t33;
	t34 = sin(qJ(5));
	t35 = cos(qJ(5));
	t22 = t32 * t35 + t34 * t38;
	t20 = 0.1e1 / t22 ^ 2;
	t21 = t32 * t34 - t35 * t38;
	t36 = t21 ^ 2 * t20 + 0.1e1;
	t24 = cos(t25);
	t19 = 0.1e1 / t36;
	t18 = (0.1e1 + t39) * t32 / (t32 ^ 2 * t39 + 0.1e1);
	t17 = -t23 * t37 + t24 * t29;
	t16 = 0.1e1 / t17 ^ 2;
	t1 = [0, 0, t18, 0, 0; 0, 0, (-t29 / t17 - (t32 * t40 + t24 * t30 + (-t24 * t37 - t40) * t18) * t30 * t16) * t33 / (t33 ^ 2 * t28 * t16 + 0.1e1), 0, 0; 0, 0, (-t35 / t22 - t34 * t21 * t20) * t33 * t30 * t19, 0, t36 * t19;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end