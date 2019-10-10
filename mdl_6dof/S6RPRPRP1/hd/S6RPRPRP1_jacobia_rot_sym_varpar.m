% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP1
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
%   Wie in S6RPRPRP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:29
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:29:04
	% EndTime: 2019-10-10 00:29:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:29:04
	% EndTime: 2019-10-10 00:29:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:29:04
	% EndTime: 2019-10-10 00:29:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:29:04
	% EndTime: 2019-10-10 00:29:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:29:04
	% EndTime: 2019-10-10 00:29:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:29:04
	% EndTime: 2019-10-10 00:29:04
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (347->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t49 = qJ(3) + pkin(10);
	t47 = cos(t49);
	t45 = sin(t49);
	t50 = qJ(1) + pkin(9);
	t46 = sin(t50);
	t58 = t46 * t45;
	t40 = atan2(-t58, -t47);
	t38 = sin(t40);
	t39 = cos(t40);
	t31 = -t38 * t58 - t39 * t47;
	t30 = 0.1e1 / t31 ^ 2;
	t48 = cos(t50);
	t64 = t30 * t48 ^ 2;
	t52 = cos(qJ(5));
	t54 = t48 * t52;
	t51 = sin(qJ(5));
	t57 = t46 * t51;
	t37 = t47 * t54 + t57;
	t35 = 0.1e1 / t37 ^ 2;
	t55 = t48 * t51;
	t56 = t46 * t52;
	t36 = t47 * t55 - t56;
	t63 = t35 * t36;
	t62 = t38 * t47;
	t42 = t45 ^ 2;
	t61 = t42 / t47 ^ 2;
	t60 = t45 * t48;
	t41 = 0.1e1 / (t46 ^ 2 * t61 + 0.1e1);
	t59 = t46 * t41;
	t53 = t36 ^ 2 * t35 + 0.1e1;
	t43 = 0.1e1 / t47;
	t34 = 0.1e1 / t37;
	t33 = 0.1e1 / t53;
	t32 = (0.1e1 + t61) * t59;
	t29 = 0.1e1 / t31;
	t28 = 0.1e1 / (t42 * t64 + 0.1e1);
	t1 = [t43 * t41 * t60, 0, t32, 0, 0, 0; (-t29 * t58 - (-t39 * t42 * t43 * t59 + (t41 - 0.1e1) * t45 * t38) * t45 * t64) * t28, 0, (t47 * t29 - (-t46 * t62 + t39 * t45 + (-t39 * t58 + t62) * t32) * t45 * t30) * t48 * t28, 0, 0, 0; ((-t47 * t57 - t54) * t34 - (-t47 * t56 + t55) * t63) * t33, 0, (-t34 * t51 + t52 * t63) * t33 * t60, 0, t53 * t33, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:29:04
	% EndTime: 2019-10-10 00:29:04
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (347->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t53 = qJ(3) + pkin(10);
	t51 = cos(t53);
	t49 = sin(t53);
	t54 = qJ(1) + pkin(9);
	t50 = sin(t54);
	t62 = t50 * t49;
	t44 = atan2(-t62, -t51);
	t42 = sin(t44);
	t43 = cos(t44);
	t35 = -t42 * t62 - t43 * t51;
	t34 = 0.1e1 / t35 ^ 2;
	t52 = cos(t54);
	t68 = t34 * t52 ^ 2;
	t56 = cos(qJ(5));
	t58 = t52 * t56;
	t55 = sin(qJ(5));
	t61 = t50 * t55;
	t41 = t51 * t58 + t61;
	t39 = 0.1e1 / t41 ^ 2;
	t59 = t52 * t55;
	t60 = t50 * t56;
	t40 = t51 * t59 - t60;
	t67 = t39 * t40;
	t66 = t42 * t51;
	t46 = t49 ^ 2;
	t65 = t46 / t51 ^ 2;
	t64 = t49 * t52;
	t45 = 0.1e1 / (t50 ^ 2 * t65 + 0.1e1);
	t63 = t50 * t45;
	t57 = t40 ^ 2 * t39 + 0.1e1;
	t47 = 0.1e1 / t51;
	t38 = 0.1e1 / t41;
	t37 = 0.1e1 / t57;
	t36 = (0.1e1 + t65) * t63;
	t33 = 0.1e1 / t35;
	t32 = 0.1e1 / (t46 * t68 + 0.1e1);
	t1 = [t47 * t45 * t64, 0, t36, 0, 0, 0; (-t33 * t62 - (-t43 * t46 * t47 * t63 + (t45 - 0.1e1) * t49 * t42) * t49 * t68) * t32, 0, (t51 * t33 - (-t50 * t66 + t43 * t49 + (-t43 * t62 + t66) * t36) * t49 * t34) * t52 * t32, 0, 0, 0; ((-t51 * t61 - t58) * t38 - (-t51 * t60 + t59) * t67) * t37, 0, (-t38 * t55 + t56 * t67) * t37 * t64, 0, t57 * t37, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end