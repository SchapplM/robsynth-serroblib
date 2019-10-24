% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRPR1
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
%   Wie in S5PPRPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:18
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PPRPR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:05
	% EndTime: 2019-10-24 10:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:05
	% EndTime: 2019-10-24 10:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:05
	% EndTime: 2019-10-24 10:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:05
	% EndTime: 2019-10-24 10:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:05
	% EndTime: 2019-10-24 10:18:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (120->15), mult. (108->32), div. (27->8), fcn. (157->9), ass. (0->22)
	t27 = pkin(8) + qJ(3);
	t26 = cos(t27);
	t25 = sin(t27);
	t29 = sin(pkin(7));
	t32 = t29 * t25;
	t21 = atan2(-t32, -t26);
	t19 = sin(t21);
	t35 = t19 * t26;
	t23 = t25 ^ 2;
	t34 = t23 / t26 ^ 2;
	t31 = cos(pkin(7));
	t33 = t26 * t31;
	t30 = cos(pkin(9));
	t28 = sin(pkin(9));
	t20 = cos(t21);
	t18 = t29 * t28 + t30 * t33;
	t17 = t28 * t33 - t29 * t30;
	t16 = 0.1e1 / t18 ^ 2;
	t14 = (0.1e1 + t34) * t29 / (t29 ^ 2 * t34 + 0.1e1);
	t13 = -t19 * t32 - t20 * t26;
	t12 = 0.1e1 / t13 ^ 2;
	t1 = [0, 0, t14, 0, 0; 0, 0, (t26 / t13 - (-t29 * t35 + t20 * t25 + (-t20 * t32 + t35) * t14) * t25 * t12) * t31 / (t31 ^ 2 * t23 * t12 + 0.1e1), 0, 0; 0, 0, (-t28 / t18 + t30 * t17 * t16) * t31 * t25 / (t17 ^ 2 * t16 + 0.1e1), 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:05
	% EndTime: 2019-10-24 10:18:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (167->16), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->25)
	t37 = pkin(8) + qJ(3);
	t35 = cos(t37);
	t33 = sin(t37);
	t38 = sin(pkin(7));
	t41 = t38 * t33;
	t28 = atan2(-t41, -t35);
	t26 = sin(t28);
	t44 = t26 * t35;
	t30 = t33 ^ 2;
	t43 = t30 / t35 ^ 2;
	t39 = cos(pkin(7));
	t42 = t35 * t39;
	t36 = pkin(9) + qJ(5);
	t32 = sin(t36);
	t34 = cos(t36);
	t25 = t38 * t32 + t34 * t42;
	t23 = 0.1e1 / t25 ^ 2;
	t24 = t32 * t42 - t38 * t34;
	t40 = t24 ^ 2 * t23 + 0.1e1;
	t27 = cos(t28);
	t22 = (0.1e1 + t43) * t38 / (t38 ^ 2 * t43 + 0.1e1);
	t21 = 0.1e1 / t40;
	t20 = -t26 * t41 - t27 * t35;
	t19 = 0.1e1 / t20 ^ 2;
	t1 = [0, 0, t22, 0, 0; 0, 0, (t35 / t20 - (-t38 * t44 + t27 * t33 + (-t27 * t41 + t44) * t22) * t33 * t19) * t39 / (t39 ^ 2 * t30 * t19 + 0.1e1), 0, 0; 0, 0, (-t32 / t25 + t34 * t24 * t23) * t39 * t33 * t21, 0, t40 * t21;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end