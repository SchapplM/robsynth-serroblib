% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPPR2
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
%   Wie in S5PRPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:22
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRPPR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:30
	% EndTime: 2019-10-24 10:22:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:30
	% EndTime: 2019-10-24 10:22:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:30
	% EndTime: 2019-10-24 10:22:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:30
	% EndTime: 2019-10-24 10:22:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:30
	% EndTime: 2019-10-24 10:22:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (120->15), mult. (108->32), div. (27->8), fcn. (157->9), ass. (0->22)
	t28 = qJ(2) + pkin(8);
	t27 = cos(t28);
	t26 = sin(t28);
	t30 = sin(pkin(7));
	t33 = t30 * t26;
	t22 = atan2(-t33, -t27);
	t20 = sin(t22);
	t36 = t20 * t27;
	t24 = t26 ^ 2;
	t35 = t24 / t27 ^ 2;
	t32 = cos(pkin(7));
	t34 = t27 * t32;
	t31 = cos(pkin(9));
	t29 = sin(pkin(9));
	t21 = cos(t22);
	t19 = t30 * t29 + t31 * t34;
	t18 = t29 * t34 - t30 * t31;
	t17 = 0.1e1 / t19 ^ 2;
	t15 = (0.1e1 + t35) * t30 / (t30 ^ 2 * t35 + 0.1e1);
	t14 = -t20 * t33 - t21 * t27;
	t13 = 0.1e1 / t14 ^ 2;
	t1 = [0, t15, 0, 0, 0; 0, (t27 / t14 - (-t30 * t36 + t21 * t26 + (-t21 * t33 + t36) * t15) * t26 * t13) * t32 / (t32 ^ 2 * t24 * t13 + 0.1e1), 0, 0, 0; 0, (-t29 / t19 + t31 * t18 * t17) * t32 * t26 / (t18 ^ 2 * t17 + 0.1e1), 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:30
	% EndTime: 2019-10-24 10:22:30
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (167->16), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->25)
	t36 = qJ(2) + pkin(8);
	t34 = cos(t36);
	t32 = sin(t36);
	t37 = sin(pkin(7));
	t40 = t37 * t32;
	t27 = atan2(-t40, -t34);
	t25 = sin(t27);
	t43 = t25 * t34;
	t29 = t32 ^ 2;
	t42 = t29 / t34 ^ 2;
	t38 = cos(pkin(7));
	t41 = t34 * t38;
	t35 = pkin(9) + qJ(5);
	t31 = sin(t35);
	t33 = cos(t35);
	t24 = t37 * t31 + t33 * t41;
	t22 = 0.1e1 / t24 ^ 2;
	t23 = t31 * t41 - t37 * t33;
	t39 = t23 ^ 2 * t22 + 0.1e1;
	t26 = cos(t27);
	t21 = (0.1e1 + t42) * t37 / (t37 ^ 2 * t42 + 0.1e1);
	t20 = 0.1e1 / t39;
	t19 = -t25 * t40 - t26 * t34;
	t18 = 0.1e1 / t19 ^ 2;
	t1 = [0, t21, 0, 0, 0; 0, (t34 / t19 - (-t37 * t43 + t26 * t32 + (-t26 * t40 + t43) * t21) * t32 * t18) * t38 / (t38 ^ 2 * t29 * t18 + 0.1e1), 0, 0, 0; 0, (-t31 / t24 + t33 * t23 * t22) * t38 * t32 * t20, 0, 0, t39 * t20;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end