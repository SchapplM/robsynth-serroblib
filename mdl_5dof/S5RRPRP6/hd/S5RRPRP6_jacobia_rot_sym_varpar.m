% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRPRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:47
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRPRP6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:47:50
	% EndTime: 2019-12-29 18:47:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:47:50
	% EndTime: 2019-12-29 18:47:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:47:50
	% EndTime: 2019-12-29 18:47:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:47:50
	% EndTime: 2019-12-29 18:47:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:47:50
	% EndTime: 2019-12-29 18:47:50
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (238->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
	t46 = qJ(2) + pkin(8);
	t45 = cos(t46);
	t44 = sin(t46);
	t48 = sin(qJ(1));
	t56 = t48 * t44;
	t39 = atan2(-t56, -t45);
	t35 = sin(t39);
	t36 = cos(t39);
	t30 = -t35 * t56 - t36 * t45;
	t29 = 0.1e1 / t30 ^ 2;
	t50 = cos(qJ(1));
	t62 = t29 * t50 ^ 2;
	t49 = cos(qJ(4));
	t52 = t50 * t49;
	t47 = sin(qJ(4));
	t55 = t48 * t47;
	t38 = t45 * t52 + t55;
	t34 = 0.1e1 / t38 ^ 2;
	t53 = t50 * t47;
	t54 = t48 * t49;
	t37 = t45 * t53 - t54;
	t61 = t34 * t37;
	t60 = t35 * t45;
	t41 = t44 ^ 2;
	t59 = t41 / t45 ^ 2;
	t58 = t44 * t50;
	t40 = 0.1e1 / (t48 ^ 2 * t59 + 0.1e1);
	t57 = t48 * t40;
	t51 = t37 ^ 2 * t34 + 0.1e1;
	t42 = 0.1e1 / t45;
	t33 = 0.1e1 / t38;
	t32 = 0.1e1 / t51;
	t31 = (0.1e1 + t59) * t57;
	t28 = 0.1e1 / t30;
	t27 = 0.1e1 / (t41 * t62 + 0.1e1);
	t1 = [t42 * t40 * t58, t31, 0, 0, 0; (-t28 * t56 - (-t36 * t41 * t42 * t57 + (t40 - 0.1e1) * t44 * t35) * t44 * t62) * t27, (t45 * t28 - (-t48 * t60 + t36 * t44 + (-t36 * t56 + t60) * t31) * t44 * t29) * t50 * t27, 0, 0, 0; ((-t45 * t55 - t52) * t33 - (-t45 * t54 + t53) * t61) * t32, (-t33 * t47 + t49 * t61) * t32 * t58, 0, t51 * t32, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:47:43
	% EndTime: 2019-12-29 18:47:43
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (238->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
	t50 = qJ(2) + pkin(8);
	t49 = cos(t50);
	t48 = sin(t50);
	t52 = sin(qJ(1));
	t60 = t52 * t48;
	t43 = atan2(-t60, -t49);
	t39 = sin(t43);
	t40 = cos(t43);
	t34 = -t39 * t60 - t40 * t49;
	t33 = 0.1e1 / t34 ^ 2;
	t54 = cos(qJ(1));
	t66 = t33 * t54 ^ 2;
	t53 = cos(qJ(4));
	t56 = t54 * t53;
	t51 = sin(qJ(4));
	t59 = t52 * t51;
	t42 = t49 * t56 + t59;
	t38 = 0.1e1 / t42 ^ 2;
	t57 = t54 * t51;
	t58 = t52 * t53;
	t41 = t49 * t57 - t58;
	t65 = t38 * t41;
	t64 = t39 * t49;
	t45 = t48 ^ 2;
	t63 = t45 / t49 ^ 2;
	t62 = t48 * t54;
	t44 = 0.1e1 / (t52 ^ 2 * t63 + 0.1e1);
	t61 = t52 * t44;
	t55 = t41 ^ 2 * t38 + 0.1e1;
	t46 = 0.1e1 / t49;
	t37 = 0.1e1 / t42;
	t36 = 0.1e1 / t55;
	t35 = (0.1e1 + t63) * t61;
	t32 = 0.1e1 / t34;
	t31 = 0.1e1 / (t45 * t66 + 0.1e1);
	t1 = [t46 * t44 * t62, t35, 0, 0, 0; (-t32 * t60 - (-t40 * t45 * t46 * t61 + (t44 - 0.1e1) * t48 * t39) * t48 * t66) * t31, (t49 * t32 - (-t52 * t64 + t40 * t48 + (-t40 * t60 + t64) * t35) * t48 * t33) * t54 * t31, 0, 0, 0; ((-t49 * t59 - t56) * t37 - (-t49 * t58 + t57) * t65) * t36, (-t37 * t51 + t53 * t65) * t36 * t62, 0, t55 * t36, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end