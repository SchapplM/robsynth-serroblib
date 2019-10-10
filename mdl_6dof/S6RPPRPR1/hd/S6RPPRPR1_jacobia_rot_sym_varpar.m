% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR1
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
%   Wie in S6RPPRPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:34
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (316->22), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->36)
	t38 = pkin(10) + qJ(4);
	t36 = cos(t38);
	t34 = sin(t38);
	t39 = qJ(1) + pkin(9);
	t35 = sin(t39);
	t46 = t35 * t34;
	t29 = atan2(-t46, -t36);
	t27 = sin(t29);
	t28 = cos(t29);
	t20 = -t27 * t46 - t28 * t36;
	t19 = 0.1e1 / t20 ^ 2;
	t37 = cos(t39);
	t52 = t19 * t37 ^ 2;
	t41 = cos(pkin(11));
	t42 = t37 * t41;
	t40 = sin(pkin(11));
	t45 = t35 * t40;
	t26 = t36 * t42 + t45;
	t24 = 0.1e1 / t26 ^ 2;
	t43 = t37 * t40;
	t44 = t35 * t41;
	t25 = t36 * t43 - t44;
	t51 = t24 * t25;
	t50 = t27 * t36;
	t31 = t34 ^ 2;
	t49 = t31 / t36 ^ 2;
	t48 = t34 * t37;
	t30 = 0.1e1 / (t35 ^ 2 * t49 + 0.1e1);
	t47 = t35 * t30;
	t32 = 0.1e1 / t36;
	t23 = 0.1e1 / t26;
	t22 = 0.1e1 / (t25 ^ 2 * t24 + 0.1e1);
	t21 = (0.1e1 + t49) * t47;
	t18 = 0.1e1 / t20;
	t17 = 0.1e1 / (t31 * t52 + 0.1e1);
	t1 = [t32 * t30 * t48, 0, 0, t21, 0, 0; (-t18 * t46 - (-t28 * t31 * t32 * t47 + (t30 - 0.1e1) * t34 * t27) * t34 * t52) * t17, 0, 0, (t36 * t18 - (-t35 * t50 + t28 * t34 + (-t28 * t46 + t50) * t21) * t34 * t19) * t37 * t17, 0, 0; ((-t36 * t45 - t42) * t23 - (-t36 * t44 + t43) * t51) * t22, 0, 0, (-t23 * t40 + t41 * t51) * t22 * t48, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:34:17
	% EndTime: 2019-10-09 23:34:17
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (395->23), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->36)
	t52 = pkin(10) + qJ(4);
	t49 = cos(t52);
	t46 = sin(t52);
	t53 = qJ(1) + pkin(9);
	t47 = sin(t53);
	t58 = t47 * t46;
	t40 = atan2(-t58, -t49);
	t38 = sin(t40);
	t39 = cos(t40);
	t31 = -t38 * t58 - t39 * t49;
	t30 = 0.1e1 / t31 ^ 2;
	t50 = cos(t53);
	t63 = t30 * t50 ^ 2;
	t51 = pkin(11) + qJ(6);
	t45 = sin(t51);
	t48 = cos(t51);
	t55 = t50 * t48;
	t37 = t47 * t45 + t49 * t55;
	t35 = 0.1e1 / t37 ^ 2;
	t56 = t50 * t45;
	t36 = -t47 * t48 + t49 * t56;
	t62 = t35 * t36;
	t42 = t46 ^ 2;
	t61 = t42 / t49 ^ 2;
	t60 = t46 * t50;
	t41 = 0.1e1 / (t47 ^ 2 * t61 + 0.1e1);
	t59 = t47 * t41;
	t57 = t47 * t49;
	t54 = t36 ^ 2 * t35 + 0.1e1;
	t43 = 0.1e1 / t49;
	t34 = 0.1e1 / t37;
	t33 = (0.1e1 + t61) * t59;
	t32 = 0.1e1 / t54;
	t29 = 0.1e1 / t31;
	t28 = 0.1e1 / (t42 * t63 + 0.1e1);
	t1 = [t43 * t41 * t60, 0, 0, t33, 0, 0; (-t29 * t58 - (-t39 * t42 * t43 * t59 + (t41 - 0.1e1) * t46 * t38) * t46 * t63) * t28, 0, 0, (t49 * t29 - (-t38 * t57 + t39 * t46 + (t38 * t49 - t39 * t58) * t33) * t46 * t30) * t50 * t28, 0, 0; ((-t45 * t57 - t55) * t34 - (-t48 * t57 + t56) * t62) * t32, 0, 0, (-t34 * t45 + t48 * t62) * t32 * t60, 0, t54 * t32;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end