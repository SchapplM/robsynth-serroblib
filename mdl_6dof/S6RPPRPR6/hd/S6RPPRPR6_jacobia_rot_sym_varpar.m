% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RPPRPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:42
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (36->14), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->27)
	t30 = cos(qJ(1));
	t26 = t30 ^ 2;
	t27 = sin(qJ(4));
	t28 = sin(qJ(1));
	t29 = cos(qJ(4));
	t34 = t28 * t29;
	t21 = atan2(t34, t27);
	t17 = sin(t21);
	t18 = cos(t21);
	t15 = t17 * t34 + t18 * t27;
	t14 = 0.1e1 / t15 ^ 2;
	t39 = t14 * t29;
	t38 = t17 * t27;
	t25 = t29 ^ 2;
	t31 = t27 ^ 2;
	t37 = 0.1e1 / t31 * t25;
	t32 = t28 ^ 2;
	t36 = 0.1e1 / t32 * t26;
	t19 = 0.1e1 / (t32 * t37 + 0.1e1);
	t35 = t28 * t19;
	t33 = t29 * t30;
	t22 = 0.1e1 / t27;
	t20 = 0.1e1 / (t31 * t36 + 0.1e1);
	t16 = (-0.1e1 - t37) * t35;
	t13 = 0.1e1 / t15;
	t12 = 0.1e1 / (t26 * t25 * t14 + 0.1e1);
	t1 = [t22 * t19 * t33, 0, 0, t16, 0, 0; (t13 * t34 + (t18 * t22 * t25 * t35 + (-t19 + 0.1e1) * t29 * t17) * t26 * t39) * t12, 0, 0, (t27 * t13 + (-t28 * t38 + t18 * t29 + (t18 * t34 - t38) * t16) * t39) * t30 * t12, 0, 0; (0.1e1 + t36) * t27 * t20, 0, 0, -0.1e1 / t28 * t20 * t33, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:42:44
	% EndTime: 2019-10-09 23:42:44
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (87->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->33)
	t43 = cos(qJ(4));
	t40 = sin(qJ(4));
	t41 = sin(qJ(1));
	t49 = t41 * t40;
	t33 = atan2(-t49, t43);
	t31 = sin(t33);
	t32 = cos(t33);
	t23 = -t31 * t49 + t32 * t43;
	t22 = 0.1e1 / t23 ^ 2;
	t44 = cos(qJ(1));
	t54 = t22 * t44 ^ 2;
	t42 = cos(qJ(6));
	t39 = sin(qJ(6));
	t47 = t44 * t39;
	t30 = -t41 * t42 - t43 * t47;
	t27 = 0.1e1 / t30 ^ 2;
	t46 = t44 * t42;
	t28 = t41 * t39 - t43 * t46;
	t53 = t27 * t28;
	t36 = t40 ^ 2;
	t52 = t36 / t43 ^ 2;
	t51 = t40 * t44;
	t34 = 0.1e1 / (t41 ^ 2 * t52 + 0.1e1);
	t50 = t41 * t34;
	t48 = t41 * t43;
	t45 = t28 ^ 2 * t27 + 0.1e1;
	t37 = 0.1e1 / t43;
	t26 = 0.1e1 / t30;
	t25 = (-0.1e1 - t52) * t50;
	t24 = 0.1e1 / t45;
	t21 = 0.1e1 / t23;
	t20 = 0.1e1 / (t36 * t54 + 0.1e1);
	t1 = [-t37 * t34 * t51, 0, 0, t25, 0, 0; (-t21 * t49 - (t32 * t36 * t37 * t50 + (t34 - 0.1e1) * t40 * t31) * t40 * t54) * t20, 0, 0, (t43 * t21 - (-t31 * t48 - t32 * t40 + (-t31 * t43 - t32 * t49) * t25) * t40 * t22) * t44 * t20, 0, 0; ((-t42 * t48 - t47) * t26 + (t39 * t48 - t46) * t53) * t24, 0, 0, (-t26 * t42 + t39 * t53) * t24 * t51, 0, t45 * t24;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end