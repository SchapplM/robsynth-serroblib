% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR4
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
%   Wie in S6RPRPPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:20
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (221->21), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->35)
	t35 = pkin(9) + qJ(3);
	t34 = cos(t35);
	t33 = sin(t35);
	t38 = sin(qJ(1));
	t44 = t38 * t33;
	t28 = atan2(-t44, -t34);
	t26 = sin(t28);
	t27 = cos(t28);
	t19 = -t26 * t44 - t27 * t34;
	t18 = 0.1e1 / t19 ^ 2;
	t39 = cos(qJ(1));
	t50 = t18 * t39 ^ 2;
	t37 = cos(pkin(10));
	t40 = t39 * t37;
	t36 = sin(pkin(10));
	t43 = t38 * t36;
	t25 = t34 * t40 + t43;
	t23 = 0.1e1 / t25 ^ 2;
	t41 = t39 * t36;
	t42 = t38 * t37;
	t24 = t34 * t41 - t42;
	t49 = t23 * t24;
	t48 = t26 * t34;
	t30 = t33 ^ 2;
	t47 = t30 / t34 ^ 2;
	t46 = t33 * t39;
	t29 = 0.1e1 / (t38 ^ 2 * t47 + 0.1e1);
	t45 = t38 * t29;
	t31 = 0.1e1 / t34;
	t22 = 0.1e1 / t25;
	t21 = 0.1e1 / (t24 ^ 2 * t23 + 0.1e1);
	t20 = (0.1e1 + t47) * t45;
	t17 = 0.1e1 / t19;
	t16 = 0.1e1 / (t30 * t50 + 0.1e1);
	t1 = [t31 * t29 * t46, 0, t20, 0, 0, 0; (-t17 * t44 - (-t27 * t30 * t31 * t45 + (t29 - 0.1e1) * t33 * t26) * t33 * t50) * t16, 0, (t34 * t17 - (-t38 * t48 + t27 * t33 + (-t27 * t44 + t48) * t20) * t33 * t18) * t39 * t16, 0, 0, 0; ((-t34 * t43 - t40) * t22 - (-t34 * t42 + t41) * t49) * t21, 0, (-t22 * t36 + t37 * t49) * t21 * t46, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (227->21), mult. (329->51), div. (61->11), fcn. (494->9), ass. (0->36)
	t49 = pkin(9) + qJ(3);
	t46 = sin(t49);
	t66 = t46 ^ 2;
	t47 = cos(t49);
	t52 = cos(pkin(10));
	t54 = cos(qJ(1));
	t56 = t54 * t52;
	t51 = sin(pkin(10));
	t53 = sin(qJ(1));
	t59 = t53 * t51;
	t38 = t47 * t59 + t56;
	t60 = t46 * t51;
	t34 = atan2(-t38, t60);
	t31 = sin(t34);
	t32 = cos(t34);
	t30 = -t31 * t38 + t32 * t60;
	t29 = 0.1e1 / t30 ^ 2;
	t57 = t54 * t51;
	t58 = t53 * t52;
	t40 = t47 * t57 - t58;
	t65 = t29 * t40;
	t63 = t32 * t38;
	t62 = t40 ^ 2 * t29;
	t48 = 0.1e1 / t51;
	t61 = 0.1e1 / t46 * t48;
	t41 = t47 * t56 + t59;
	t37 = 0.1e1 / t41 ^ 2;
	t55 = t54 ^ 2 * t66 * t37;
	t45 = 0.1e1 / t66;
	t36 = 0.1e1 / t41;
	t35 = 0.1e1 / (0.1e1 + t55);
	t33 = 0.1e1 / (0.1e1 + t38 ^ 2 * t45 / t51 ^ 2);
	t28 = 0.1e1 / t30;
	t27 = (t38 * t45 * t47 * t48 + t53) * t33;
	t26 = 0.1e1 / (0.1e1 + t62);
	t1 = [-t40 * t33 * t61, 0, t27, 0, 0, 0; (-t38 * t28 - (-t31 + (t61 * t63 + t31) * t33) * t62) * t26, 0, (t27 * t63 * t65 + (-t54 * t46 * t28 - (t32 * t47 + (-t27 + t53) * t46 * t31) * t65) * t51) * t26, 0, 0, 0; (t53 * t36 + (-t47 * t58 + t57) * t54 * t37) * t46 * t35, 0, (-t36 * t47 * t54 - t52 * t55) * t35, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:20:28
	% EndTime: 2019-10-10 00:20:28
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (265->25), mult. (347->66), div. (52->9), fcn. (503->11), ass. (0->42)
	t55 = pkin(9) + qJ(3);
	t54 = cos(t55);
	t53 = sin(t55);
	t59 = sin(qJ(1));
	t67 = t59 * t53;
	t49 = atan2(t67, t54);
	t46 = sin(t49);
	t47 = cos(t49);
	t36 = t46 * t67 + t47 * t54;
	t35 = 0.1e1 / t36 ^ 2;
	t61 = cos(qJ(1));
	t73 = t35 * t61 ^ 2;
	t56 = sin(pkin(10));
	t64 = t61 * t56;
	t57 = cos(pkin(10));
	t65 = t59 * t57;
	t44 = t54 * t64 - t65;
	t63 = t61 * t57;
	t66 = t59 * t56;
	t45 = t54 * t63 + t66;
	t58 = sin(qJ(6));
	t60 = cos(qJ(6));
	t41 = t44 * t58 + t45 * t60;
	t39 = 0.1e1 / t41 ^ 2;
	t40 = -t44 * t60 + t45 * t58;
	t72 = t39 * t40;
	t71 = t46 * t54;
	t50 = t53 ^ 2;
	t70 = t50 / t54 ^ 2;
	t69 = t53 * t61;
	t48 = 0.1e1 / (t59 ^ 2 * t70 + 0.1e1);
	t68 = t59 * t48;
	t62 = t40 ^ 2 * t39 + 0.1e1;
	t51 = 0.1e1 / t54;
	t43 = -t54 * t65 + t64;
	t42 = -t54 * t66 - t63;
	t38 = 0.1e1 / t41;
	t37 = (0.1e1 + t70) * t68;
	t34 = 0.1e1 / t36;
	t33 = 0.1e1 / (t50 * t73 + 0.1e1);
	t32 = 0.1e1 / t62;
	t1 = [t51 * t48 * t69, 0, t37, 0, 0, 0; (t34 * t67 + (t47 * t50 * t51 * t68 + (-t48 + 0.1e1) * t53 * t46) * t53 * t73) * t33, 0, (-t54 * t34 + (t59 * t71 - t47 * t53 + (t47 * t67 - t71) * t37) * t53 * t35) * t61 * t33, 0, 0, 0; ((-t42 * t60 + t43 * t58) * t38 - (t42 * t58 + t43 * t60) * t72) * t32, 0, ((t56 * t60 - t57 * t58) * t38 - (-t56 * t58 - t57 * t60) * t72) * t32 * t69, 0, 0, t62 * t32;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end