% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPRP1
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
%   Wie in S5PRPRP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRPRP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:29:18
	% EndTime: 2019-12-05 15:29:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:29:18
	% EndTime: 2019-12-05 15:29:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:29:18
	% EndTime: 2019-12-05 15:29:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:29:18
	% EndTime: 2019-12-05 15:29:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:29:18
	% EndTime: 2019-12-05 15:29:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:29:18
	% EndTime: 2019-12-05 15:29:18
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (264->18), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->28)
	t35 = pkin(7) + qJ(2);
	t33 = cos(t35);
	t36 = t33 ^ 2;
	t34 = pkin(8) + qJ(4);
	t32 = cos(t34);
	t30 = sin(t34);
	t31 = sin(t35);
	t37 = t31 * t30;
	t22 = atan2(-t37, -t32);
	t20 = sin(t22);
	t21 = cos(t22);
	t18 = -t20 * t37 - t21 * t32;
	t17 = 0.1e1 / t18 ^ 2;
	t42 = t17 * t30;
	t41 = t20 * t32;
	t25 = t30 ^ 2;
	t28 = 0.1e1 / t32 ^ 2;
	t40 = t25 * t28;
	t26 = t31 ^ 2;
	t39 = t26 / t36;
	t23 = 0.1e1 / (t26 * t40 + 0.1e1);
	t38 = t31 * t23;
	t27 = 0.1e1 / t32;
	t24 = 0.1e1 / (t28 * t39 + 0.1e1);
	t19 = (0.1e1 + t40) * t38;
	t16 = 0.1e1 / t18;
	t15 = 0.1e1 / (t36 * t25 * t17 + 0.1e1);
	t1 = [0, t33 * t30 * t27 * t23, 0, t19, 0; 0, (-t16 * t37 - (-t21 * t25 * t27 * t38 + (t23 - 0.1e1) * t30 * t20) * t36 * t42) * t15, 0, (t32 * t16 - (-t31 * t41 + t21 * t30 + (-t21 * t37 + t41) * t19) * t42) * t33 * t15, 0; 0, (-0.1e1 - t39) * t27 * t24, 0, -0.1e1 / t33 * t28 * t24 * t37, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end