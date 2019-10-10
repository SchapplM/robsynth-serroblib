% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP1
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
%   Wie in S6RPRRPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:09
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
	% DurationCPUTime: 0.15s
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
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:42
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (266->22), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
	t47 = qJ(1) + pkin(9);
	t45 = cos(t47);
	t62 = t45 ^ 2;
	t52 = cos(qJ(3));
	t43 = sin(t47);
	t51 = sin(qJ(3));
	t58 = t43 * t51;
	t40 = atan2(-t58, -t52);
	t38 = sin(t40);
	t39 = cos(t40);
	t32 = -t38 * t58 - t39 * t52;
	t31 = 0.1e1 / t32 ^ 2;
	t61 = t31 * t51;
	t46 = qJ(4) + pkin(10);
	t42 = sin(t46);
	t44 = cos(t46);
	t55 = t45 * t52;
	t37 = t43 * t42 + t44 * t55;
	t35 = 0.1e1 / t37 ^ 2;
	t36 = t42 * t55 - t43 * t44;
	t60 = t35 * t36;
	t48 = t51 ^ 2;
	t54 = t48 / t52 ^ 2;
	t41 = 0.1e1 / (t43 ^ 2 * t54 + 0.1e1);
	t59 = t43 * t41;
	t57 = t43 * t52;
	t56 = t45 * t51;
	t53 = t36 ^ 2 * t35 + 0.1e1;
	t49 = 0.1e1 / t52;
	t34 = 0.1e1 / t37;
	t33 = (0.1e1 + t54) * t59;
	t30 = 0.1e1 / t32;
	t29 = 0.1e1 / t53;
	t28 = 0.1e1 / (t62 * t48 * t31 + 0.1e1);
	t1 = [t49 * t41 * t56, 0, t33, 0, 0, 0; (-t30 * t58 - (-t39 * t48 * t49 * t59 + (t41 - 0.1e1) * t51 * t38) * t62 * t61) * t28, 0, (t52 * t30 - (-t38 * t57 + t39 * t51 + (t38 * t52 - t39 * t58) * t33) * t61) * t45 * t28, 0, 0, 0; ((-t42 * t57 - t45 * t44) * t34 - (t45 * t42 - t44 * t57) * t60) * t29, 0, (-t34 * t42 + t44 * t60) * t29 * t56, t53 * t29, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:42
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (665->28), mult. (509->70), div. (100->11), fcn. (770->9), ass. (0->40)
	t65 = sin(qJ(3));
	t77 = t65 ^ 2;
	t60 = qJ(4) + pkin(10);
	t56 = sin(t60);
	t58 = cos(t60);
	t61 = qJ(1) + pkin(9);
	t59 = cos(t61);
	t57 = sin(t61);
	t66 = cos(qJ(3));
	t71 = t57 * t66;
	t46 = t56 * t71 + t59 * t58;
	t68 = t65 * t56;
	t43 = atan2(-t46, t68);
	t39 = sin(t43);
	t40 = cos(t43);
	t38 = -t39 * t46 + t40 * t68;
	t37 = 0.1e1 / t38 ^ 2;
	t69 = t59 * t66;
	t49 = t56 * t69 - t57 * t58;
	t76 = t37 * t49;
	t74 = t40 * t46;
	t73 = t49 ^ 2 * t37;
	t53 = 0.1e1 / t56;
	t63 = 0.1e1 / t65;
	t72 = t53 * t63;
	t70 = t59 * t65;
	t50 = t57 * t56 + t58 * t69;
	t45 = 0.1e1 / t50 ^ 2;
	t67 = t59 ^ 2 * t77 * t45;
	t64 = 0.1e1 / t77;
	t54 = 0.1e1 / t56 ^ 2;
	t48 = -t59 * t56 + t58 * t71;
	t44 = 0.1e1 / t50;
	t42 = 0.1e1 / (t46 ^ 2 * t64 * t54 + 0.1e1);
	t41 = 0.1e1 / (0.1e1 + t67);
	t36 = 0.1e1 / t38;
	t35 = (t46 * t53 * t64 * t66 + t57) * t42;
	t34 = 0.1e1 / (0.1e1 + t73);
	t33 = (t46 * t54 * t58 - t48 * t53) * t63 * t42;
	t1 = [-t49 * t42 * t72, 0, t35, t33, 0, 0; (-t46 * t36 - (-t39 + (t72 * t74 + t39) * t42) * t73) * t34, 0, (t35 * t74 * t76 + (-t36 * t70 - (t40 * t66 + (-t35 + t57) * t65 * t39) * t76) * t56) * t34, (t50 * t36 - (t40 * t65 * t58 - t39 * t48 + (-t39 * t68 - t74) * t33) * t76) * t34, 0, 0; (-t45 * t48 * t59 + t44 * t57) * t65 * t41, 0, (-t44 * t69 - t58 * t67) * t41, -t49 * t45 * t41 * t70, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end