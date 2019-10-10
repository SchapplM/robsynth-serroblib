% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP3
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
%   Wie in S6RPRRPP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:13
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:10
	% EndTime: 2019-10-10 01:13:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:10
	% EndTime: 2019-10-10 01:13:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:10
	% EndTime: 2019-10-10 01:13:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:10
	% EndTime: 2019-10-10 01:13:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:10
	% EndTime: 2019-10-10 01:13:10
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (218->21), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
	t36 = qJ(1) + pkin(9);
	t35 = cos(t36);
	t54 = t35 ^ 2;
	t43 = cos(qJ(3));
	t34 = sin(t36);
	t41 = sin(qJ(3));
	t49 = t34 * t41;
	t32 = atan2(-t49, -t43);
	t30 = sin(t32);
	t31 = cos(t32);
	t23 = -t30 * t49 - t31 * t43;
	t22 = 0.1e1 / t23 ^ 2;
	t53 = t22 * t41;
	t40 = sin(qJ(4));
	t42 = cos(qJ(4));
	t45 = t42 * t43;
	t29 = t34 * t40 + t35 * t45;
	t27 = 0.1e1 / t29 ^ 2;
	t46 = t40 * t43;
	t28 = -t34 * t42 + t35 * t46;
	t52 = t27 * t28;
	t51 = t30 * t43;
	t37 = t41 ^ 2;
	t47 = t37 / t43 ^ 2;
	t33 = 0.1e1 / (t34 ^ 2 * t47 + 0.1e1);
	t50 = t34 * t33;
	t48 = t35 * t41;
	t44 = t28 ^ 2 * t27 + 0.1e1;
	t38 = 0.1e1 / t43;
	t26 = 0.1e1 / t29;
	t25 = (0.1e1 + t47) * t50;
	t24 = 0.1e1 / t44;
	t21 = 0.1e1 / t23;
	t20 = 0.1e1 / (t54 * t37 * t22 + 0.1e1);
	t1 = [t38 * t33 * t48, 0, t25, 0, 0, 0; (-t21 * t49 - (-t31 * t37 * t38 * t50 + (t33 - 0.1e1) * t41 * t30) * t54 * t53) * t20, 0, (t43 * t21 - (-t34 * t51 + t31 * t41 + (-t31 * t49 + t51) * t25) * t53) * t35 * t20, 0, 0, 0; ((-t34 * t46 - t35 * t42) * t26 - (-t34 * t45 + t35 * t40) * t52) * t24, 0, (-t26 * t40 + t42 * t52) * t24 * t48, t44 * t24, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:10
	% EndTime: 2019-10-10 01:13:10
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (369->27), mult. (486->71), div. (108->11), fcn. (755->9), ass. (0->38)
	t46 = qJ(1) + pkin(9);
	t44 = sin(t46);
	t45 = cos(t46);
	t53 = cos(qJ(4));
	t51 = sin(qJ(4));
	t54 = cos(qJ(3));
	t57 = t51 * t54;
	t35 = t44 * t57 + t45 * t53;
	t52 = sin(qJ(3));
	t56 = t52 * t51;
	t34 = atan2(-t35, t56);
	t31 = sin(t34);
	t32 = cos(t34);
	t29 = -t31 * t35 + t32 * t56;
	t28 = 0.1e1 / t29 ^ 2;
	t38 = -t44 * t53 + t45 * t57;
	t64 = t28 * t38;
	t55 = t53 * t54;
	t39 = t44 * t51 + t45 * t55;
	t43 = 0.1e1 / t45 ^ 2;
	t50 = 0.1e1 / t52 ^ 2;
	t30 = 0.1e1 / (t39 ^ 2 * t43 * t50 + 0.1e1);
	t49 = 0.1e1 / t52;
	t63 = t30 * t49;
	t61 = t32 * t35;
	t60 = t38 ^ 2 * t28;
	t47 = 0.1e1 / t51;
	t59 = t47 * t49;
	t58 = t50 * t54;
	t48 = 0.1e1 / t51 ^ 2;
	t42 = 0.1e1 / t45;
	t37 = t44 * t55 - t45 * t51;
	t33 = 0.1e1 / (t35 ^ 2 * t50 * t48 + 0.1e1);
	t27 = 0.1e1 / t29;
	t26 = (t35 * t47 * t58 + t44) * t33;
	t25 = 0.1e1 / (0.1e1 + t60);
	t24 = (t35 * t48 * t53 - t37 * t47) * t49 * t33;
	t1 = [-t38 * t33 * t59, 0, t26, t24, 0, 0; (-t35 * t27 - (-t31 + (t59 * t61 + t31) * t33) * t60) * t25, 0, (t26 * t61 * t64 + (-t45 * t52 * t27 - (t32 * t54 + (-t26 + t44) * t52 * t31) * t64) * t51) * t25, (t39 * t27 - (t32 * t52 * t53 - t31 * t37 + (-t31 * t56 - t61) * t24) * t64) * t25, 0, 0; (t39 * t43 * t44 - t37 * t42) * t63, 0, (-t39 * t42 * t58 - t53) * t30, -t38 * t42 * t63, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:13:10
	% EndTime: 2019-10-10 01:13:10
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (369->27), mult. (486->71), div. (108->11), fcn. (755->9), ass. (0->38)
	t45 = qJ(1) + pkin(9);
	t43 = sin(t45);
	t44 = cos(t45);
	t50 = sin(qJ(4));
	t52 = cos(qJ(4));
	t53 = cos(qJ(3));
	t54 = t52 * t53;
	t36 = t43 * t54 - t44 * t50;
	t51 = sin(qJ(3));
	t55 = t51 * t52;
	t34 = atan2(-t36, t55);
	t31 = sin(t34);
	t32 = cos(t34);
	t29 = -t31 * t36 + t32 * t55;
	t28 = 0.1e1 / t29 ^ 2;
	t39 = t43 * t50 + t44 * t54;
	t63 = t28 * t39;
	t56 = t50 * t53;
	t38 = t43 * t52 - t44 * t56;
	t42 = 0.1e1 / t44 ^ 2;
	t47 = 0.1e1 / t51 ^ 2;
	t30 = 0.1e1 / (t38 ^ 2 * t42 * t47 + 0.1e1);
	t46 = 0.1e1 / t51;
	t62 = t30 * t46;
	t60 = t32 * t36;
	t59 = t39 ^ 2 * t28;
	t48 = 0.1e1 / t52;
	t58 = t46 * t48;
	t57 = t47 * t53;
	t49 = 0.1e1 / t52 ^ 2;
	t41 = 0.1e1 / t44;
	t35 = t43 * t56 + t44 * t52;
	t33 = 0.1e1 / (t36 ^ 2 * t47 * t49 + 0.1e1);
	t27 = 0.1e1 / t29;
	t26 = (t36 * t48 * t57 + t43) * t33;
	t25 = 0.1e1 / (0.1e1 + t59);
	t24 = (-t36 * t49 * t50 + t35 * t48) * t46 * t33;
	t1 = [-t39 * t33 * t58, 0, t26, t24, 0, 0; (-t36 * t27 - (-t31 + (t58 * t60 + t31) * t33) * t59) * t25, 0, (t26 * t60 * t63 + (-t44 * t51 * t27 - (t32 * t53 + (-t26 + t43) * t31 * t51) * t63) * t52) * t25, (t38 * t27 - (-t32 * t51 * t50 + t31 * t35 + (-t31 * t55 - t60) * t24) * t63) * t25, 0, 0; (t38 * t42 * t43 + t35 * t41) * t62, 0, (-t38 * t41 * t57 + t50) * t30, -t39 * t41 * t62, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end