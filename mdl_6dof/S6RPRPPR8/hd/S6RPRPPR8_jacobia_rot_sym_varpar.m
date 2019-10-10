% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR8
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
%   Wie in S6RPRPPR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (58->14), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->26)
	t27 = sin(qJ(1));
	t30 = t27 ^ 2;
	t26 = sin(qJ(3));
	t28 = cos(qJ(3));
	t29 = cos(qJ(1));
	t31 = t29 * t28;
	t18 = atan2(-t31, t26);
	t16 = sin(t18);
	t17 = cos(t18);
	t14 = -t16 * t31 + t17 * t26;
	t13 = 0.1e1 / t14 ^ 2;
	t36 = t13 * t28;
	t35 = t16 * t26;
	t22 = 0.1e1 / t26 ^ 2;
	t24 = t28 ^ 2;
	t34 = t22 * t24;
	t25 = t29 ^ 2;
	t33 = 0.1e1 / t30 * t25;
	t20 = 0.1e1 / (t25 * t34 + 0.1e1);
	t32 = t29 * t20;
	t21 = 0.1e1 / t26;
	t19 = 0.1e1 / (t22 * t33 + 0.1e1);
	t15 = (0.1e1 + t34) * t32;
	t12 = 0.1e1 / t14;
	t11 = 0.1e1 / (t30 * t24 * t13 + 0.1e1);
	t1 = [t27 * t28 * t21 * t20, 0, t15, 0, 0, 0; (-t12 * t31 + (-t17 * t21 * t24 * t32 + (-t20 + 0.1e1) * t28 * t16) * t30 * t36) * t11, 0, (t26 * t12 + (t29 * t35 + t17 * t28 + (-t17 * t31 - t35) * t15) * t36) * t27 * t11, 0, 0, 0; (0.1e1 + t33) * t21 * t19, 0, 0.1e1 / t27 * t22 * t19 * t31, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:27:23
	% EndTime: 2019-10-10 00:27:23
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (64->19), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->34)
	t43 = cos(qJ(3));
	t40 = sin(qJ(3));
	t44 = cos(qJ(1));
	t47 = t44 * t40;
	t34 = atan2(t47, t43);
	t31 = sin(t34);
	t32 = cos(t34);
	t23 = t31 * t47 + t32 * t43;
	t22 = 0.1e1 / t23 ^ 2;
	t41 = sin(qJ(1));
	t55 = t22 * t41 ^ 2;
	t42 = cos(qJ(6));
	t39 = sin(qJ(6));
	t48 = t44 * t39;
	t50 = t41 * t43;
	t30 = -t42 * t50 - t48;
	t27 = 0.1e1 / t30 ^ 2;
	t46 = t44 * t42;
	t28 = t39 * t50 - t46;
	t54 = t27 * t28;
	t53 = t31 * t43;
	t36 = t40 ^ 2;
	t52 = t36 / t43 ^ 2;
	t51 = t40 * t41;
	t33 = 0.1e1 / (t44 ^ 2 * t52 + 0.1e1);
	t49 = t44 * t33;
	t45 = t28 ^ 2 * t27 + 0.1e1;
	t37 = 0.1e1 / t43;
	t26 = 0.1e1 / t30;
	t25 = (0.1e1 + t52) * t49;
	t24 = 0.1e1 / t45;
	t21 = 0.1e1 / t23;
	t20 = 0.1e1 / (t36 * t55 + 0.1e1);
	t1 = [-t37 * t33 * t51, 0, t25, 0, 0, 0; (t21 * t47 - (-t32 * t36 * t37 * t49 + (t33 - 0.1e1) * t40 * t31) * t40 * t55) * t20, 0, (t43 * t21 - (t44 * t53 - t32 * t40 + (t32 * t47 - t53) * t25) * t40 * t22) * t41 * t20, 0, 0, 0; ((-t41 * t42 - t43 * t48) * t26 + (t41 * t39 - t43 * t46) * t54) * t24, 0, (t26 * t39 + t42 * t54) * t24 * t51, 0, 0, t45 * t24;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end