% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP8
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
%   Wie in S6RRPRRP8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:41
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (99->20), mult. (197->54), div. (47->9), fcn. (297->9), ass. (0->32)
	t34 = cos(qJ(2));
	t32 = sin(qJ(2));
	t33 = sin(qJ(1));
	t39 = t33 * t32;
	t25 = atan2(-t39, -t34);
	t23 = sin(t25);
	t24 = cos(t25);
	t16 = -t23 * t39 - t24 * t34;
	t15 = 0.1e1 / t16 ^ 2;
	t35 = cos(qJ(1));
	t44 = t15 * t35 ^ 2;
	t30 = sin(pkin(10));
	t31 = cos(pkin(10));
	t36 = t35 * t31;
	t22 = t33 * t30 + t34 * t36;
	t20 = 0.1e1 / t22 ^ 2;
	t37 = t35 * t30;
	t21 = -t33 * t31 + t34 * t37;
	t43 = t20 * t21;
	t27 = t32 ^ 2;
	t42 = t27 / t34 ^ 2;
	t41 = t32 * t35;
	t26 = 0.1e1 / (t33 ^ 2 * t42 + 0.1e1);
	t40 = t33 * t26;
	t38 = t33 * t34;
	t28 = 0.1e1 / t34;
	t19 = 0.1e1 / t22;
	t18 = (0.1e1 + t42) * t40;
	t17 = 0.1e1 / (t21 ^ 2 * t20 + 0.1e1);
	t14 = 0.1e1 / t16;
	t13 = 0.1e1 / (t27 * t44 + 0.1e1);
	t1 = [t28 * t26 * t41, t18, 0, 0, 0, 0; (-t14 * t39 - (-t24 * t27 * t28 * t40 + (t26 - 0.1e1) * t32 * t23) * t32 * t44) * t13, (t34 * t14 - (-t23 * t38 + t24 * t32 + (t23 * t34 - t24 * t39) * t18) * t32 * t15) * t35 * t13, 0, 0, 0, 0; ((-t30 * t38 - t36) * t19 - (-t31 * t38 + t37) * t43) * t17, (-t19 * t30 + t31 * t43) * t17 * t41, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (157->21), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->34)
	t46 = cos(qJ(2));
	t44 = sin(qJ(2));
	t45 = sin(qJ(1));
	t52 = t45 * t44;
	t36 = atan2(-t52, -t46);
	t34 = sin(t36);
	t35 = cos(t36);
	t28 = -t34 * t52 - t35 * t46;
	t27 = 0.1e1 / t28 ^ 2;
	t47 = cos(qJ(1));
	t57 = t27 * t47 ^ 2;
	t40 = pkin(10) + qJ(4);
	t38 = sin(t40);
	t39 = cos(t40);
	t49 = t47 * t39;
	t33 = t45 * t38 + t46 * t49;
	t31 = 0.1e1 / t33 ^ 2;
	t50 = t47 * t38;
	t32 = -t45 * t39 + t46 * t50;
	t56 = t31 * t32;
	t41 = t44 ^ 2;
	t55 = t41 / t46 ^ 2;
	t54 = t44 * t47;
	t37 = 0.1e1 / (t45 ^ 2 * t55 + 0.1e1);
	t53 = t45 * t37;
	t51 = t45 * t46;
	t48 = t32 ^ 2 * t31 + 0.1e1;
	t42 = 0.1e1 / t46;
	t30 = 0.1e1 / t33;
	t29 = (0.1e1 + t55) * t53;
	t26 = 0.1e1 / t28;
	t25 = 0.1e1 / t48;
	t24 = 0.1e1 / (t41 * t57 + 0.1e1);
	t1 = [t42 * t37 * t54, t29, 0, 0, 0, 0; (-t26 * t52 - (-t35 * t41 * t42 * t53 + (t37 - 0.1e1) * t44 * t34) * t44 * t57) * t24, (t46 * t26 - (-t34 * t51 + t35 * t44 + (t34 * t46 - t35 * t52) * t29) * t44 * t27) * t47 * t24, 0, 0, 0, 0; ((-t38 * t51 - t49) * t30 - (-t39 * t51 + t50) * t56) * t25, (-t30 * t38 + t39 * t56) * t25 * t54, 0, t48 * t25, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (243->21), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->35)
	t59 = cos(qJ(2));
	t57 = sin(qJ(2));
	t58 = sin(qJ(1));
	t65 = t58 * t57;
	t49 = atan2(-t65, -t59);
	t47 = sin(t49);
	t48 = cos(t49);
	t41 = -t47 * t65 - t48 * t59;
	t40 = 0.1e1 / t41 ^ 2;
	t60 = cos(qJ(1));
	t70 = t40 * t60 ^ 2;
	t53 = pkin(10) + qJ(4) + qJ(5);
	t51 = sin(t53);
	t52 = cos(t53);
	t62 = t60 * t52;
	t46 = t58 * t51 + t59 * t62;
	t44 = 0.1e1 / t46 ^ 2;
	t63 = t60 * t51;
	t45 = -t58 * t52 + t59 * t63;
	t69 = t44 * t45;
	t54 = t57 ^ 2;
	t68 = t54 / t59 ^ 2;
	t67 = t57 * t60;
	t50 = 0.1e1 / (t58 ^ 2 * t68 + 0.1e1);
	t66 = t58 * t50;
	t64 = t58 * t59;
	t61 = t45 ^ 2 * t44 + 0.1e1;
	t55 = 0.1e1 / t59;
	t43 = 0.1e1 / t46;
	t42 = (0.1e1 + t68) * t66;
	t39 = 0.1e1 / t41;
	t38 = 0.1e1 / (t54 * t70 + 0.1e1);
	t37 = 0.1e1 / t61;
	t36 = t61 * t37;
	t1 = [t55 * t50 * t67, t42, 0, 0, 0, 0; (-t39 * t65 - (-t48 * t54 * t55 * t66 + (t50 - 0.1e1) * t57 * t47) * t57 * t70) * t38, (t59 * t39 - (-t47 * t64 + t48 * t57 + (t47 * t59 - t48 * t65) * t42) * t57 * t40) * t60 * t38, 0, 0, 0, 0; ((-t51 * t64 - t62) * t43 - (-t52 * t64 + t63) * t69) * t37, (-t43 * t51 + t52 * t69) * t37 * t67, 0, t36, t36, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:04
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (1004->28), mult. (689->68), div. (139->11), fcn. (1046->9), ass. (0->42)
	t77 = sin(qJ(2));
	t92 = t77 ^ 2;
	t72 = pkin(10) + qJ(4) + qJ(5);
	t70 = sin(t72);
	t71 = cos(t72);
	t80 = cos(qJ(1));
	t82 = t80 * t71;
	t78 = sin(qJ(1));
	t79 = cos(qJ(2));
	t84 = t78 * t79;
	t61 = t70 * t84 + t82;
	t86 = t77 * t70;
	t57 = atan2(-t61, t86);
	t54 = sin(t57);
	t55 = cos(t57);
	t52 = -t54 * t61 + t55 * t86;
	t51 = 0.1e1 / t52 ^ 2;
	t83 = t80 * t70;
	t64 = -t78 * t71 + t79 * t83;
	t91 = t51 * t64;
	t89 = t55 * t61;
	t88 = t64 ^ 2 * t51;
	t68 = 0.1e1 / t70;
	t74 = 0.1e1 / t77;
	t87 = t68 * t74;
	t85 = t77 * t80;
	t65 = t78 * t70 + t79 * t82;
	t60 = 0.1e1 / t65 ^ 2;
	t81 = t80 ^ 2 * t92 * t60;
	t75 = 0.1e1 / t92;
	t69 = 0.1e1 / t70 ^ 2;
	t63 = t71 * t84 - t83;
	t59 = 0.1e1 / t65;
	t58 = 0.1e1 / (0.1e1 + t81);
	t56 = 0.1e1 / (t61 ^ 2 * t75 * t69 + 0.1e1);
	t53 = t64 * t60 * t58 * t85;
	t50 = 0.1e1 / t52;
	t49 = (t61 * t68 * t75 * t79 + t78) * t56;
	t48 = 0.1e1 / (0.1e1 + t88);
	t47 = (t61 * t69 * t71 - t63 * t68) * t74 * t56;
	t46 = (t65 * t50 - (t55 * t77 * t71 - t54 * t63 + (-t54 * t86 - t89) * t47) * t91) * t48;
	t1 = [-t64 * t56 * t87, t49, 0, t47, t47, 0; (-t61 * t50 - (-t54 + (t87 * t89 + t54) * t56) * t88) * t48, (t49 * t89 * t91 + (-t50 * t85 - (t55 * t79 + (-t49 + t78) * t77 * t54) * t91) * t70) * t48, 0, t46, t46, 0; (-t60 * t63 * t80 + t59 * t78) * t77 * t58, (-t59 * t79 * t80 - t71 * t81) * t58, 0, -t53, -t53, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end