% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR8
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
%   Wie in S6RPRPRR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:21
	% EndTime: 2019-10-10 00:58:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:21
	% EndTime: 2019-10-10 00:58:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:21
	% EndTime: 2019-10-10 00:58:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:21
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (215->20), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
	t46 = qJ(3) + pkin(10);
	t44 = sin(t46);
	t45 = cos(t46);
	t50 = cos(qJ(1));
	t54 = t50 * t45;
	t39 = atan2(-t54, t44);
	t37 = sin(t39);
	t38 = cos(t39);
	t30 = -t37 * t54 + t38 * t44;
	t29 = 0.1e1 / t30 ^ 2;
	t48 = sin(qJ(1));
	t62 = t29 * t48 ^ 2;
	t47 = sin(qJ(5));
	t53 = t50 * t47;
	t49 = cos(qJ(5));
	t56 = t48 * t49;
	t36 = t44 * t56 + t53;
	t34 = 0.1e1 / t36 ^ 2;
	t52 = t50 * t49;
	t57 = t48 * t47;
	t35 = t44 * t57 - t52;
	t61 = t34 * t35;
	t60 = t37 * t44;
	t43 = t45 ^ 2;
	t59 = 0.1e1 / t44 ^ 2 * t43;
	t58 = t45 * t48;
	t40 = 0.1e1 / (t50 ^ 2 * t59 + 0.1e1);
	t55 = t50 * t40;
	t51 = t35 ^ 2 * t34 + 0.1e1;
	t41 = 0.1e1 / t44;
	t33 = 0.1e1 / t36;
	t32 = 0.1e1 / t51;
	t31 = (0.1e1 + t59) * t55;
	t28 = 0.1e1 / t30;
	t27 = 0.1e1 / (t43 * t62 + 0.1e1);
	t1 = [t41 * t40 * t58, 0, t31, 0, 0, 0; (-t28 * t54 + (-t38 * t41 * t43 * t55 + (-t40 + 0.1e1) * t45 * t37) * t45 * t62) * t27, 0, (t44 * t28 + (t50 * t60 + t38 * t45 + (-t38 * t54 - t60) * t31) * t45 * t29) * t48 * t27, 0, 0, 0; ((t44 * t53 + t56) * t33 - (t44 * t52 - t57) * t61) * t32, 0, (t33 * t47 - t49 * t61) * t32 * t58, 0, t51 * t32, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:58:22
	% EndTime: 2019-10-10 00:58:22
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (294->21), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->38)
	t63 = qJ(3) + pkin(10);
	t59 = sin(t63);
	t60 = cos(t63);
	t66 = cos(qJ(1));
	t70 = t66 * t60;
	t54 = atan2(-t70, t59);
	t52 = sin(t54);
	t53 = cos(t54);
	t45 = -t52 * t70 + t53 * t59;
	t44 = 0.1e1 / t45 ^ 2;
	t65 = sin(qJ(1));
	t78 = t44 * t65 ^ 2;
	t64 = qJ(5) + qJ(6);
	t61 = sin(t64);
	t69 = t66 * t61;
	t62 = cos(t64);
	t72 = t65 * t62;
	t51 = t59 * t72 + t69;
	t49 = 0.1e1 / t51 ^ 2;
	t68 = t66 * t62;
	t73 = t65 * t61;
	t50 = t59 * t73 - t68;
	t77 = t49 * t50;
	t76 = t52 * t59;
	t58 = t60 ^ 2;
	t75 = 0.1e1 / t59 ^ 2 * t58;
	t74 = t60 * t65;
	t55 = 0.1e1 / (t66 ^ 2 * t75 + 0.1e1);
	t71 = t66 * t55;
	t67 = t50 ^ 2 * t49 + 0.1e1;
	t56 = 0.1e1 / t59;
	t48 = 0.1e1 / t51;
	t47 = 0.1e1 / t67;
	t46 = (0.1e1 + t75) * t71;
	t43 = 0.1e1 / t45;
	t42 = 0.1e1 / (t58 * t78 + 0.1e1);
	t41 = t67 * t47;
	t1 = [t56 * t55 * t74, 0, t46, 0, 0, 0; (-t43 * t70 + (-t53 * t56 * t58 * t71 + (-t55 + 0.1e1) * t60 * t52) * t60 * t78) * t42, 0, (t59 * t43 + (t66 * t76 + t53 * t60 + (-t53 * t70 - t76) * t46) * t60 * t44) * t65 * t42, 0, 0, 0; ((t59 * t69 + t72) * t48 - (t59 * t68 - t73) * t77) * t47, 0, (t48 * t61 - t62 * t77) * t47 * t74, 0, t41, t41;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end