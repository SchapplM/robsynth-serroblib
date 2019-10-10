% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP8
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
%   Wie in S6RPRRPP8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:21
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:48
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:48
	% EndTime: 2019-10-10 01:21:48
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:48
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (86->19), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->36)
	t40 = sin(qJ(1));
	t56 = t40 ^ 2;
	t39 = sin(qJ(3));
	t42 = cos(qJ(3));
	t43 = cos(qJ(1));
	t45 = t43 * t42;
	t33 = atan2(-t45, t39);
	t31 = sin(t33);
	t32 = cos(t33);
	t24 = -t31 * t45 + t32 * t39;
	t23 = 0.1e1 / t24 ^ 2;
	t55 = t23 * t42;
	t38 = sin(qJ(4));
	t47 = t43 * t38;
	t41 = cos(qJ(4));
	t50 = t40 * t41;
	t30 = t39 * t50 + t47;
	t28 = 0.1e1 / t30 ^ 2;
	t46 = t43 * t41;
	t51 = t40 * t38;
	t29 = t39 * t51 - t46;
	t54 = t28 * t29;
	t53 = t31 * t39;
	t37 = t42 ^ 2;
	t52 = 0.1e1 / t39 ^ 2 * t37;
	t49 = t40 * t42;
	t34 = 0.1e1 / (t43 ^ 2 * t52 + 0.1e1);
	t48 = t43 * t34;
	t44 = t29 ^ 2 * t28 + 0.1e1;
	t35 = 0.1e1 / t39;
	t27 = 0.1e1 / t30;
	t26 = (0.1e1 + t52) * t48;
	t25 = 0.1e1 / t44;
	t22 = 0.1e1 / t24;
	t21 = 0.1e1 / (t56 * t37 * t23 + 0.1e1);
	t1 = [t35 * t34 * t49, 0, t26, 0, 0, 0; (-t22 * t45 + (-t32 * t35 * t37 * t48 + (-t34 + 0.1e1) * t42 * t31) * t56 * t55) * t21, 0, (t39 * t22 + (t43 * t53 + t32 * t42 + (-t32 * t45 - t53) * t26) * t55) * t40 * t21, 0, 0, 0; ((t39 * t47 + t50) * t27 - (t39 * t46 - t51) * t54) * t25, 0, (t27 * t38 - t41 * t54) * t25 * t49, t44 * t25, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (157->24), mult. (486->69), div. (108->11), fcn. (755->9), ass. (0->39)
	t50 = sin(qJ(3));
	t49 = sin(qJ(4));
	t54 = cos(qJ(1));
	t56 = t54 * t49;
	t51 = sin(qJ(1));
	t52 = cos(qJ(4));
	t58 = t51 * t52;
	t40 = t50 * t56 + t58;
	t53 = cos(qJ(3));
	t57 = t53 * t49;
	t37 = atan2(t40, t57);
	t33 = sin(t37);
	t34 = cos(t37);
	t32 = t33 * t40 + t34 * t57;
	t31 = 0.1e1 / t32 ^ 2;
	t55 = t54 * t52;
	t59 = t51 * t49;
	t38 = t50 * t59 - t55;
	t66 = t31 * t38;
	t64 = t34 * t40;
	t39 = t50 * t58 + t56;
	t46 = 0.1e1 / t51 ^ 2;
	t48 = 0.1e1 / t53 ^ 2;
	t35 = 0.1e1 / (t39 ^ 2 * t46 * t48 + 0.1e1);
	t47 = 0.1e1 / t53;
	t63 = t35 * t47;
	t62 = t38 ^ 2 * t31;
	t43 = 0.1e1 / t49;
	t61 = t43 * t47;
	t60 = t48 * t50;
	t45 = 0.1e1 / t51;
	t44 = 0.1e1 / t49 ^ 2;
	t41 = t50 * t55 - t59;
	t36 = 0.1e1 / (t40 ^ 2 * t48 * t44 + 0.1e1);
	t30 = 0.1e1 / t32;
	t29 = (t40 * t43 * t60 + t54) * t36;
	t28 = 0.1e1 / (0.1e1 + t62);
	t27 = (-t40 * t44 * t52 + t41 * t43) * t47 * t36;
	t1 = [-t38 * t36 * t61, 0, t29, t27, 0, 0; (t40 * t30 - (-t33 + (-t61 * t64 + t33) * t36) * t62) * t28, 0, (-t29 * t64 * t66 + (t51 * t53 * t30 - (-t34 * t50 + (-t29 + t54) * t53 * t33) * t66) * t49) * t28, (t39 * t30 - (t34 * t53 * t52 + t33 * t41 + (-t33 * t57 + t64) * t27) * t66) * t28, 0, 0; (t39 * t46 * t54 - t41 * t45) * t63, 0, (-t39 * t45 * t60 - t52) * t35, t38 * t45 * t63, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:21:49
	% EndTime: 2019-10-10 01:21:49
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (157->24), mult. (486->69), div. (108->11), fcn. (755->9), ass. (0->39)
	t50 = sin(qJ(3));
	t52 = cos(qJ(4));
	t54 = cos(qJ(1));
	t55 = t54 * t52;
	t49 = sin(qJ(4));
	t51 = sin(qJ(1));
	t59 = t51 * t49;
	t41 = t50 * t55 - t59;
	t53 = cos(qJ(3));
	t57 = t53 * t52;
	t37 = atan2(t41, t57);
	t33 = sin(t37);
	t34 = cos(t37);
	t32 = t33 * t41 + t34 * t57;
	t31 = 0.1e1 / t32 ^ 2;
	t56 = t54 * t49;
	t58 = t51 * t52;
	t39 = t50 * t58 + t56;
	t66 = t31 * t39;
	t64 = t34 * t41;
	t38 = -t50 * t59 + t55;
	t44 = 0.1e1 / t51 ^ 2;
	t48 = 0.1e1 / t53 ^ 2;
	t35 = 0.1e1 / (t38 ^ 2 * t44 * t48 + 0.1e1);
	t47 = 0.1e1 / t53;
	t63 = t35 * t47;
	t62 = t39 ^ 2 * t31;
	t45 = 0.1e1 / t52;
	t61 = t45 * t47;
	t60 = t48 * t50;
	t46 = 0.1e1 / t52 ^ 2;
	t43 = 0.1e1 / t51;
	t40 = t50 * t56 + t58;
	t36 = 0.1e1 / (t41 ^ 2 * t48 * t46 + 0.1e1);
	t30 = 0.1e1 / t32;
	t29 = (t41 * t45 * t60 + t54) * t36;
	t28 = 0.1e1 / (0.1e1 + t62);
	t27 = (t41 * t46 * t49 - t40 * t45) * t47 * t36;
	t1 = [-t39 * t36 * t61, 0, t29, t27, 0, 0; (t41 * t30 - (-t33 + (-t61 * t64 + t33) * t36) * t62) * t28, 0, (-t29 * t64 * t66 + (t51 * t53 * t30 - (-t34 * t50 + (-t29 + t54) * t53 * t33) * t66) * t52) * t28, (t38 * t30 - (-t34 * t53 * t49 - t33 * t40 + (-t33 * t57 + t64) * t27) * t66) * t28, 0, 0; (t38 * t44 * t54 + t40 * t43) * t63, 0, (-t38 * t43 * t60 + t49) * t35, t39 * t43 * t63, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end