% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR1
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
%   Wie in S6RPRPPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:15
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:14
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (316->22), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->36)
	t41 = qJ(3) + pkin(10);
	t39 = cos(t41);
	t37 = sin(t41);
	t42 = qJ(1) + pkin(9);
	t38 = sin(t42);
	t49 = t38 * t37;
	t32 = atan2(-t49, -t39);
	t30 = sin(t32);
	t31 = cos(t32);
	t23 = -t30 * t49 - t31 * t39;
	t22 = 0.1e1 / t23 ^ 2;
	t40 = cos(t42);
	t55 = t22 * t40 ^ 2;
	t44 = cos(pkin(11));
	t45 = t40 * t44;
	t43 = sin(pkin(11));
	t48 = t38 * t43;
	t29 = t39 * t45 + t48;
	t27 = 0.1e1 / t29 ^ 2;
	t46 = t40 * t43;
	t47 = t38 * t44;
	t28 = t39 * t46 - t47;
	t54 = t27 * t28;
	t53 = t30 * t39;
	t34 = t37 ^ 2;
	t52 = t34 / t39 ^ 2;
	t51 = t37 * t40;
	t33 = 0.1e1 / (t38 ^ 2 * t52 + 0.1e1);
	t50 = t38 * t33;
	t35 = 0.1e1 / t39;
	t26 = 0.1e1 / t29;
	t25 = 0.1e1 / (t28 ^ 2 * t27 + 0.1e1);
	t24 = (0.1e1 + t52) * t50;
	t21 = 0.1e1 / t23;
	t20 = 0.1e1 / (t34 * t55 + 0.1e1);
	t1 = [t35 * t33 * t51, 0, t24, 0, 0, 0; (-t21 * t49 - (-t31 * t34 * t35 * t50 + (t33 - 0.1e1) * t37 * t30) * t37 * t55) * t20, 0, (t39 * t21 - (-t38 * t53 + t31 * t37 + (-t31 * t49 + t53) * t24) * t37 * t22) * t40 * t20, 0, 0, 0; ((-t39 * t48 - t45) * t26 - (-t39 * t47 + t46) * t54) * t25, 0, (-t26 * t43 + t44 * t54) * t25 * t51, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:15:13
	% EndTime: 2019-10-10 00:15:13
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (395->23), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->36)
	t55 = qJ(3) + pkin(10);
	t52 = cos(t55);
	t49 = sin(t55);
	t56 = qJ(1) + pkin(9);
	t50 = sin(t56);
	t61 = t50 * t49;
	t43 = atan2(-t61, -t52);
	t41 = sin(t43);
	t42 = cos(t43);
	t34 = -t41 * t61 - t42 * t52;
	t33 = 0.1e1 / t34 ^ 2;
	t53 = cos(t56);
	t66 = t33 * t53 ^ 2;
	t54 = pkin(11) + qJ(6);
	t48 = sin(t54);
	t51 = cos(t54);
	t58 = t53 * t51;
	t40 = t50 * t48 + t52 * t58;
	t38 = 0.1e1 / t40 ^ 2;
	t59 = t53 * t48;
	t39 = -t50 * t51 + t52 * t59;
	t65 = t38 * t39;
	t45 = t49 ^ 2;
	t64 = t45 / t52 ^ 2;
	t63 = t49 * t53;
	t44 = 0.1e1 / (t50 ^ 2 * t64 + 0.1e1);
	t62 = t50 * t44;
	t60 = t50 * t52;
	t57 = t39 ^ 2 * t38 + 0.1e1;
	t46 = 0.1e1 / t52;
	t37 = 0.1e1 / t40;
	t36 = (0.1e1 + t64) * t62;
	t35 = 0.1e1 / t57;
	t32 = 0.1e1 / t34;
	t31 = 0.1e1 / (t45 * t66 + 0.1e1);
	t1 = [t46 * t44 * t63, 0, t36, 0, 0, 0; (-t32 * t61 - (-t42 * t45 * t46 * t62 + (t44 - 0.1e1) * t49 * t41) * t49 * t66) * t31, 0, (t52 * t32 - (-t41 * t60 + t42 * t49 + (t41 * t52 - t42 * t61) * t36) * t49 * t33) * t53 * t31, 0, 0, 0; ((-t48 * t60 - t58) * t37 - (-t51 * t60 + t59) * t65) * t35, 0, (-t37 * t48 + t51 * t65) * t35 * t63, 0, 0, t57 * t35;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end