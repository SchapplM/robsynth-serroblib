% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR8
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
%   Wie in S6RPPRRR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:13
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRR8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (215->20), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
	t43 = pkin(10) + qJ(4);
	t41 = sin(t43);
	t42 = cos(t43);
	t47 = cos(qJ(1));
	t51 = t47 * t42;
	t36 = atan2(-t51, t41);
	t34 = sin(t36);
	t35 = cos(t36);
	t27 = -t34 * t51 + t35 * t41;
	t26 = 0.1e1 / t27 ^ 2;
	t45 = sin(qJ(1));
	t59 = t26 * t45 ^ 2;
	t44 = sin(qJ(5));
	t50 = t47 * t44;
	t46 = cos(qJ(5));
	t53 = t45 * t46;
	t33 = t41 * t53 + t50;
	t31 = 0.1e1 / t33 ^ 2;
	t49 = t47 * t46;
	t54 = t45 * t44;
	t32 = t41 * t54 - t49;
	t58 = t31 * t32;
	t57 = t34 * t41;
	t40 = t42 ^ 2;
	t56 = 0.1e1 / t41 ^ 2 * t40;
	t55 = t42 * t45;
	t37 = 0.1e1 / (t47 ^ 2 * t56 + 0.1e1);
	t52 = t47 * t37;
	t48 = t32 ^ 2 * t31 + 0.1e1;
	t38 = 0.1e1 / t41;
	t30 = 0.1e1 / t33;
	t29 = 0.1e1 / t48;
	t28 = (0.1e1 + t56) * t52;
	t25 = 0.1e1 / t27;
	t24 = 0.1e1 / (t40 * t59 + 0.1e1);
	t1 = [t38 * t37 * t55, 0, 0, t28, 0, 0; (-t25 * t51 + (-t35 * t38 * t40 * t52 + (-t37 + 0.1e1) * t42 * t34) * t42 * t59) * t24, 0, 0, (t41 * t25 + (t47 * t57 + t35 * t42 + (-t35 * t51 - t57) * t28) * t42 * t26) * t45 * t24, 0, 0; ((t41 * t50 + t53) * t30 - (t41 * t49 - t54) * t58) * t29, 0, 0, (t30 * t44 - t46 * t58) * t29 * t55, t48 * t29, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (294->21), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->38)
	t60 = pkin(10) + qJ(4);
	t56 = sin(t60);
	t57 = cos(t60);
	t63 = cos(qJ(1));
	t67 = t63 * t57;
	t51 = atan2(-t67, t56);
	t49 = sin(t51);
	t50 = cos(t51);
	t42 = -t49 * t67 + t50 * t56;
	t41 = 0.1e1 / t42 ^ 2;
	t62 = sin(qJ(1));
	t75 = t41 * t62 ^ 2;
	t61 = qJ(5) + qJ(6);
	t58 = sin(t61);
	t66 = t63 * t58;
	t59 = cos(t61);
	t69 = t62 * t59;
	t48 = t56 * t69 + t66;
	t46 = 0.1e1 / t48 ^ 2;
	t65 = t63 * t59;
	t70 = t62 * t58;
	t47 = t56 * t70 - t65;
	t74 = t46 * t47;
	t73 = t49 * t56;
	t55 = t57 ^ 2;
	t72 = 0.1e1 / t56 ^ 2 * t55;
	t71 = t57 * t62;
	t52 = 0.1e1 / (t63 ^ 2 * t72 + 0.1e1);
	t68 = t63 * t52;
	t64 = t47 ^ 2 * t46 + 0.1e1;
	t53 = 0.1e1 / t56;
	t45 = 0.1e1 / t48;
	t44 = 0.1e1 / t64;
	t43 = (0.1e1 + t72) * t68;
	t40 = 0.1e1 / t42;
	t39 = 0.1e1 / (t55 * t75 + 0.1e1);
	t38 = t64 * t44;
	t1 = [t53 * t52 * t71, 0, 0, t43, 0, 0; (-t40 * t67 + (-t50 * t53 * t55 * t68 + (-t52 + 0.1e1) * t57 * t49) * t57 * t75) * t39, 0, 0, (t56 * t40 + (t63 * t73 + t50 * t57 + (-t50 * t67 - t73) * t43) * t57 * t41) * t62 * t39, 0, 0; ((t56 * t66 + t69) * t45 - (t56 * t65 - t70) * t74) * t44, 0, 0, (t45 * t58 - t59 * t74) * t44 * t71, t38, t38;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end