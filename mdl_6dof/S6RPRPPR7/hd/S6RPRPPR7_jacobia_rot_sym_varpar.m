% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR7
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
%   Wie in S6RPRPPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (168->15), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->28)
	t35 = sin(qJ(1));
	t33 = t35 ^ 2;
	t32 = qJ(3) + pkin(9);
	t30 = sin(t32);
	t31 = cos(t32);
	t36 = cos(qJ(1));
	t39 = t36 * t31;
	t24 = atan2(-t39, t30);
	t22 = sin(t24);
	t23 = cos(t24);
	t20 = -t22 * t39 + t23 * t30;
	t19 = 0.1e1 / t20 ^ 2;
	t45 = t19 * t31;
	t44 = t22 * t30;
	t29 = t31 ^ 2;
	t37 = t30 ^ 2;
	t43 = 0.1e1 / t37 * t29;
	t42 = t31 * t35;
	t38 = t36 ^ 2;
	t41 = t33 / t38;
	t25 = 0.1e1 / (t38 * t43 + 0.1e1);
	t40 = t36 * t25;
	t27 = 0.1e1 / t30;
	t26 = 0.1e1 / (t37 * t41 + 0.1e1);
	t21 = (0.1e1 + t43) * t40;
	t18 = 0.1e1 / t20;
	t17 = 0.1e1 / (t33 * t29 * t19 + 0.1e1);
	t1 = [t27 * t25 * t42, 0, t21, 0, 0, 0; (-t18 * t39 + (-t23 * t27 * t29 * t40 + (-t25 + 0.1e1) * t31 * t22) * t33 * t45) * t17, 0, (t30 * t18 + (t36 * t44 + t23 * t31 + (-t23 * t39 - t44) * t21) * t45) * t35 * t17, 0, 0, 0; (0.1e1 + t41) * t30 * t26, 0, 0.1e1 / t36 * t26 * t42, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:25:42
	% EndTime: 2019-10-10 00:25:42
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (193->20), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
	t46 = qJ(3) + pkin(9);
	t45 = cos(t46);
	t44 = sin(t46);
	t50 = cos(qJ(1));
	t54 = t50 * t44;
	t40 = atan2(t54, t45);
	t37 = sin(t40);
	t38 = cos(t40);
	t30 = t37 * t54 + t38 * t45;
	t29 = 0.1e1 / t30 ^ 2;
	t48 = sin(qJ(1));
	t62 = t29 * t48 ^ 2;
	t49 = cos(qJ(6));
	t52 = t50 * t49;
	t47 = sin(qJ(6));
	t57 = t48 * t47;
	t36 = -t45 * t57 + t52;
	t34 = 0.1e1 / t36 ^ 2;
	t53 = t50 * t47;
	t56 = t48 * t49;
	t35 = t45 * t56 + t53;
	t61 = t34 * t35;
	t60 = t37 * t45;
	t41 = t44 ^ 2;
	t59 = t41 / t45 ^ 2;
	t58 = t44 * t48;
	t39 = 0.1e1 / (t50 ^ 2 * t59 + 0.1e1);
	t55 = t50 * t39;
	t51 = t35 ^ 2 * t34 + 0.1e1;
	t42 = 0.1e1 / t45;
	t33 = 0.1e1 / t36;
	t32 = 0.1e1 / t51;
	t31 = (0.1e1 + t59) * t55;
	t28 = 0.1e1 / t30;
	t27 = 0.1e1 / (t41 * t62 + 0.1e1);
	t1 = [-t42 * t39 * t58, 0, t31, 0, 0, 0; (t28 * t54 - (-t38 * t41 * t42 * t55 + (t39 - 0.1e1) * t44 * t37) * t44 * t62) * t27, 0, (t45 * t28 - (t50 * t60 - t38 * t44 + (t38 * t54 - t60) * t31) * t44 * t29) * t48 * t27, 0, 0, 0; ((t45 * t52 - t57) * t33 - (-t45 * t53 - t56) * t61) * t32, 0, (-t33 * t49 - t47 * t61) * t32 * t58, 0, 0, t51 * t32;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end