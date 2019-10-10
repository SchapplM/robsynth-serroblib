% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR3
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
%   Wie in S6RPRPPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:18
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:44
	% EndTime: 2019-10-10 00:18:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:43
	% EndTime: 2019-10-10 00:18:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:43
	% EndTime: 2019-10-10 00:18:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:43
	% EndTime: 2019-10-10 00:18:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:43
	% EndTime: 2019-10-10 00:18:43
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (153->17), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->27)
	t26 = qJ(1) + pkin(9);
	t25 = cos(t26);
	t32 = t25 ^ 2;
	t31 = cos(qJ(3));
	t24 = sin(t26);
	t30 = sin(qJ(3));
	t34 = t24 * t30;
	t20 = atan2(-t34, -t31);
	t17 = sin(t20);
	t18 = cos(t20);
	t15 = -t17 * t34 - t18 * t31;
	t14 = 0.1e1 / t15 ^ 2;
	t38 = t14 * t30;
	t37 = t17 * t31;
	t22 = t24 ^ 2;
	t36 = t22 / t32;
	t27 = t30 ^ 2;
	t29 = 0.1e1 / t31 ^ 2;
	t33 = t27 * t29;
	t21 = 0.1e1 / (t22 * t33 + 0.1e1);
	t35 = t24 * t21;
	t28 = 0.1e1 / t31;
	t19 = 0.1e1 / (t29 * t36 + 0.1e1);
	t16 = (0.1e1 + t33) * t35;
	t13 = 0.1e1 / t15;
	t12 = 0.1e1 / (t32 * t27 * t14 + 0.1e1);
	t1 = [t25 * t30 * t28 * t21, 0, t16, 0, 0, 0; (-t13 * t34 - (-t18 * t27 * t28 * t35 + (t21 - 0.1e1) * t30 * t17) * t32 * t38) * t12, 0, (t31 * t13 - (-t24 * t37 + t18 * t30 + (-t18 * t34 + t37) * t16) * t38) * t25 * t12, 0, 0, 0; (-0.1e1 - t36) * t28 * t19, 0, -0.1e1 / t25 * t29 * t19 * t34, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:43
	% EndTime: 2019-10-10 00:18:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:18:43
	% EndTime: 2019-10-10 00:18:44
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (196->21), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
	t38 = qJ(1) + pkin(9);
	t37 = cos(t38);
	t56 = t37 ^ 2;
	t43 = sin(qJ(3));
	t36 = sin(t38);
	t45 = cos(qJ(3));
	t51 = t36 * t45;
	t34 = atan2(-t51, t43);
	t32 = sin(t34);
	t33 = cos(t34);
	t27 = -t32 * t51 + t33 * t43;
	t24 = 0.1e1 / t27 ^ 2;
	t55 = t24 * t45;
	t42 = sin(qJ(6));
	t44 = cos(qJ(6));
	t47 = t43 * t44;
	t31 = -t36 * t42 + t37 * t47;
	t29 = 0.1e1 / t31 ^ 2;
	t48 = t42 * t43;
	t30 = t36 * t44 + t37 * t48;
	t54 = t29 * t30;
	t53 = t32 * t43;
	t41 = t45 ^ 2;
	t49 = 0.1e1 / t43 ^ 2 * t41;
	t35 = 0.1e1 / (t36 ^ 2 * t49 + 0.1e1);
	t52 = t36 * t35;
	t50 = t37 * t45;
	t46 = t30 ^ 2 * t29 + 0.1e1;
	t39 = 0.1e1 / t43;
	t28 = 0.1e1 / t31;
	t26 = (0.1e1 + t49) * t52;
	t25 = 0.1e1 / t46;
	t23 = 0.1e1 / t27;
	t22 = 0.1e1 / (t56 * t41 * t24 + 0.1e1);
	t1 = [-t39 * t35 * t50, 0, t26, 0, 0, 0; (-t23 * t51 - (t33 * t39 * t41 * t52 + (t35 - 0.1e1) * t45 * t32) * t56 * t55) * t22, 0, (-t43 * t23 - (t36 * t53 + t33 * t45 + (-t33 * t51 - t53) * t26) * t55) * t37 * t22, 0, 0, 0; ((-t36 * t48 + t37 * t44) * t28 - (-t36 * t47 - t37 * t42) * t54) * t25, 0, (t28 * t42 - t44 * t54) * t25 * t50, 0, 0, t46 * t25;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end