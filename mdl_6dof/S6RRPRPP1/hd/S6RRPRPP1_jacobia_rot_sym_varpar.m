% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPP1
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
%   Wie in S6RRPRPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:55
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:40
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (238->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
	t46 = qJ(2) + pkin(9);
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
	t1 = [t42 * t40 * t58, t31, 0, 0, 0, 0; (-t28 * t56 - (-t36 * t41 * t42 * t57 + (t40 - 0.1e1) * t44 * t35) * t44 * t62) * t27, (t45 * t28 - (-t48 * t60 + t36 * t44 + (-t36 * t56 + t60) * t31) * t44 * t29) * t50 * t27, 0, 0, 0, 0; ((-t45 * t55 - t52) * t33 - (-t45 * t54 + t53) * t61) * t32, (-t33 * t47 + t49 * t61) * t32 * t58, 0, t51 * t32, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:40
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (286->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t54 = qJ(2) + pkin(9);
	t52 = cos(t54);
	t50 = sin(t54);
	t55 = sin(qJ(1));
	t61 = t55 * t50;
	t44 = atan2(-t61, -t52);
	t42 = sin(t44);
	t43 = cos(t44);
	t35 = -t42 * t61 - t43 * t52;
	t34 = 0.1e1 / t35 ^ 2;
	t56 = cos(qJ(1));
	t68 = t34 * t56 ^ 2;
	t53 = qJ(4) + pkin(10);
	t51 = cos(t53);
	t58 = t56 * t51;
	t49 = sin(t53);
	t62 = t55 * t49;
	t41 = t52 * t58 + t62;
	t39 = 0.1e1 / t41 ^ 2;
	t59 = t56 * t49;
	t60 = t55 * t51;
	t40 = t52 * t59 - t60;
	t67 = t39 * t40;
	t66 = t42 * t52;
	t46 = t50 ^ 2;
	t65 = t46 / t52 ^ 2;
	t64 = t50 * t56;
	t45 = 0.1e1 / (t55 ^ 2 * t65 + 0.1e1);
	t63 = t55 * t45;
	t57 = t40 ^ 2 * t39 + 0.1e1;
	t47 = 0.1e1 / t52;
	t38 = 0.1e1 / t41;
	t37 = (0.1e1 + t65) * t63;
	t36 = 0.1e1 / t57;
	t33 = 0.1e1 / t35;
	t32 = 0.1e1 / (t46 * t68 + 0.1e1);
	t1 = [t47 * t45 * t64, t37, 0, 0, 0, 0; (-t33 * t61 - (-t43 * t46 * t47 * t63 + (t45 - 0.1e1) * t50 * t42) * t50 * t68) * t32, (t52 * t33 - (-t55 * t66 + t43 * t50 + (-t43 * t61 + t66) * t37) * t50 * t34) * t56 * t32, 0, 0, 0, 0; ((-t52 * t62 - t58) * t38 - (-t52 * t60 + t59) * t67) * t36, (-t38 * t49 + t51 * t67) * t36 * t64, 0, t57 * t36, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:55:39
	% EndTime: 2019-10-10 09:55:40
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (640->28), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->42)
	t68 = qJ(2) + pkin(9);
	t64 = sin(t68);
	t84 = t64 ^ 2;
	t66 = cos(t68);
	t67 = qJ(4) + pkin(10);
	t65 = cos(t67);
	t71 = cos(qJ(1));
	t73 = t71 * t65;
	t63 = sin(t67);
	t70 = sin(qJ(1));
	t76 = t70 * t63;
	t51 = t66 * t76 + t73;
	t78 = t64 * t63;
	t47 = atan2(-t51, t78);
	t44 = sin(t47);
	t45 = cos(t47);
	t43 = -t44 * t51 + t45 * t78;
	t42 = 0.1e1 / t43 ^ 2;
	t74 = t71 * t63;
	t75 = t70 * t65;
	t54 = t66 * t74 - t75;
	t83 = t42 * t54;
	t81 = t45 * t51;
	t80 = t54 ^ 2 * t42;
	t58 = 0.1e1 / t63;
	t61 = 0.1e1 / t64;
	t79 = t58 * t61;
	t77 = t64 * t71;
	t55 = t66 * t73 + t76;
	t50 = 0.1e1 / t55 ^ 2;
	t72 = t71 ^ 2 * t84 * t50;
	t62 = 0.1e1 / t84;
	t59 = 0.1e1 / t63 ^ 2;
	t53 = t66 * t75 - t74;
	t49 = 0.1e1 / t55;
	t48 = 0.1e1 / (0.1e1 + t72);
	t46 = 0.1e1 / (t51 ^ 2 * t62 * t59 + 0.1e1);
	t41 = 0.1e1 / t43;
	t40 = (t51 * t58 * t62 * t66 + t70) * t46;
	t39 = 0.1e1 / (0.1e1 + t80);
	t38 = (t51 * t59 * t65 - t53 * t58) * t61 * t46;
	t1 = [-t54 * t46 * t79, t40, 0, t38, 0, 0; (-t51 * t41 - (-t44 + (t79 * t81 + t44) * t46) * t80) * t39, (t40 * t81 * t83 + (-t41 * t77 - (t45 * t66 + (-t40 + t70) * t64 * t44) * t83) * t63) * t39, 0, (t55 * t41 - (t45 * t64 * t65 - t44 * t53 + (-t44 * t78 - t81) * t38) * t83) * t39, 0, 0; (-t50 * t53 * t71 + t49 * t70) * t64 * t48, (-t49 * t66 * t71 - t65 * t72) * t48, 0, -t54 * t50 * t48 * t77, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end