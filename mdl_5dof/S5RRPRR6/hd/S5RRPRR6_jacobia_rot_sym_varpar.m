% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRR6
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
%   Wie in S5RRPRR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRPRR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
Ja_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (18->0), mult. (12->0), div. (10->0), fcn. (12->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (206->15), mult. (205->35), div. (45->9), fcn. (315->9), ass. (0->27)
	t51 = cos(pkin(9));
	t49 = qJ(1) + qJ(2);
	t45 = sin(t49);
	t50 = sin(pkin(9));
	t57 = t45 * t50;
	t43 = atan2(-t57, -t51);
	t41 = sin(t43);
	t42 = cos(t43);
	t35 = -t41 * t57 - t42 * t51;
	t46 = cos(t49);
	t59 = 0.1e1 / t35 ^ 2 * t46 ^ 2;
	t47 = t50 ^ 2;
	t44 = 0.1e1 / (0.1e1 + t45 ^ 2 * t47 / t51 ^ 2);
	t58 = t44 / t51;
	t52 = sin(qJ(4));
	t56 = t51 * t52;
	t53 = cos(qJ(4));
	t55 = t51 * t53;
	t40 = t45 * t52 + t46 * t55;
	t38 = 0.1e1 / t40 ^ 2;
	t39 = -t45 * t53 + t46 * t56;
	t54 = t39 ^ 2 * t38 + 0.1e1;
	t37 = t46 * t50 * t58;
	t36 = 0.1e1 / t54;
	t32 = ((-t45 * t56 - t46 * t53) / t40 - (-t45 * t55 + t46 * t52) * t39 * t38) * t36;
	t31 = (-0.1e1 / t35 * t57 - (-t42 * t45 * t47 * t58 + (t44 - 0.1e1) * t50 * t41) * t50 * t59) / (t47 * t59 + 0.1e1);
	t1 = [t37, t37, 0, 0, 0; t31, t31, 0, 0, 0; t32, t32, 0, t54 * t36, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:06
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (294->16), mult. (232->35), div. (50->9), fcn. (350->9), ass. (0->29)
	t70 = cos(pkin(9));
	t68 = qJ(1) + qJ(2);
	t62 = sin(t68);
	t69 = sin(pkin(9));
	t74 = t62 * t69;
	t59 = atan2(-t74, -t70);
	t57 = sin(t59);
	t58 = cos(t59);
	t52 = -t57 * t74 - t58 * t70;
	t64 = cos(t68);
	t76 = 0.1e1 / t52 ^ 2 * t64 ^ 2;
	t65 = t69 ^ 2;
	t60 = 0.1e1 / (0.1e1 + t62 ^ 2 * t65 / t70 ^ 2);
	t75 = t60 / t70;
	t73 = t62 * t70;
	t72 = t64 * t70;
	t67 = qJ(4) + qJ(5);
	t61 = sin(t67);
	t63 = cos(t67);
	t56 = t62 * t61 + t63 * t72;
	t53 = 0.1e1 / t56 ^ 2;
	t55 = t61 * t72 - t62 * t63;
	t71 = t55 ^ 2 * t53 + 0.1e1;
	t54 = t64 * t69 * t75;
	t50 = 0.1e1 / t71;
	t48 = t71 * t50;
	t47 = ((-t61 * t73 - t64 * t63) / t56 - (t64 * t61 - t63 * t73) * t55 * t53) * t50;
	t46 = (-0.1e1 / t52 * t74 - (-t58 * t62 * t65 * t75 + (t60 - 0.1e1) * t69 * t57) * t69 * t76) / (t65 * t76 + 0.1e1);
	t1 = [t54, t54, 0, 0, 0; t46, t46, 0, 0, 0; t47, t47, 0, t48, t48;];
	Ja_rot = t1;
end