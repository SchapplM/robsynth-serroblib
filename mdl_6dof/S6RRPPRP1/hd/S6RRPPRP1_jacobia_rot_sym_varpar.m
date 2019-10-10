% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRP1
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
%   Wie in S6RRPPRP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:27
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:26:57
	% EndTime: 2019-10-10 09:26:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:26:57
	% EndTime: 2019-10-10 09:26:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:26:57
	% EndTime: 2019-10-10 09:26:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:26:57
	% EndTime: 2019-10-10 09:26:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:26:57
	% EndTime: 2019-10-10 09:26:57
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (221->21), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->35)
	t38 = qJ(2) + pkin(9);
	t37 = cos(t38);
	t36 = sin(t38);
	t41 = sin(qJ(1));
	t47 = t41 * t36;
	t31 = atan2(-t47, -t37);
	t29 = sin(t31);
	t30 = cos(t31);
	t22 = -t29 * t47 - t30 * t37;
	t21 = 0.1e1 / t22 ^ 2;
	t42 = cos(qJ(1));
	t53 = t21 * t42 ^ 2;
	t40 = cos(pkin(10));
	t43 = t42 * t40;
	t39 = sin(pkin(10));
	t46 = t41 * t39;
	t28 = t37 * t43 + t46;
	t26 = 0.1e1 / t28 ^ 2;
	t44 = t42 * t39;
	t45 = t41 * t40;
	t27 = t37 * t44 - t45;
	t52 = t26 * t27;
	t51 = t29 * t37;
	t33 = t36 ^ 2;
	t50 = t33 / t37 ^ 2;
	t49 = t36 * t42;
	t32 = 0.1e1 / (t41 ^ 2 * t50 + 0.1e1);
	t48 = t41 * t32;
	t34 = 0.1e1 / t37;
	t25 = 0.1e1 / t28;
	t24 = 0.1e1 / (t27 ^ 2 * t26 + 0.1e1);
	t23 = (0.1e1 + t50) * t48;
	t20 = 0.1e1 / t22;
	t19 = 0.1e1 / (t33 * t53 + 0.1e1);
	t1 = [t34 * t32 * t49, t23, 0, 0, 0, 0; (-t20 * t47 - (-t30 * t33 * t34 * t48 + (t32 - 0.1e1) * t36 * t29) * t36 * t53) * t19, (t37 * t20 - (-t41 * t51 + t30 * t36 + (-t30 * t47 + t51) * t23) * t36 * t21) * t42 * t19, 0, 0, 0, 0; ((-t37 * t46 - t43) * t25 - (-t37 * t45 + t44) * t52) * t24, (-t25 * t39 + t40 * t52) * t24 * t49, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:26:57
	% EndTime: 2019-10-10 09:26:57
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (286->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t52 = qJ(2) + pkin(9);
	t50 = cos(t52);
	t48 = sin(t52);
	t53 = sin(qJ(1));
	t59 = t53 * t48;
	t42 = atan2(-t59, -t50);
	t40 = sin(t42);
	t41 = cos(t42);
	t33 = -t40 * t59 - t41 * t50;
	t32 = 0.1e1 / t33 ^ 2;
	t54 = cos(qJ(1));
	t66 = t32 * t54 ^ 2;
	t51 = pkin(10) + qJ(5);
	t49 = cos(t51);
	t56 = t54 * t49;
	t47 = sin(t51);
	t60 = t53 * t47;
	t39 = t50 * t56 + t60;
	t37 = 0.1e1 / t39 ^ 2;
	t57 = t54 * t47;
	t58 = t53 * t49;
	t38 = t50 * t57 - t58;
	t65 = t37 * t38;
	t64 = t40 * t50;
	t44 = t48 ^ 2;
	t63 = t44 / t50 ^ 2;
	t62 = t48 * t54;
	t43 = 0.1e1 / (t53 ^ 2 * t63 + 0.1e1);
	t61 = t53 * t43;
	t55 = t38 ^ 2 * t37 + 0.1e1;
	t45 = 0.1e1 / t50;
	t36 = 0.1e1 / t39;
	t35 = (0.1e1 + t63) * t61;
	t34 = 0.1e1 / t55;
	t31 = 0.1e1 / t33;
	t30 = 0.1e1 / (t44 * t66 + 0.1e1);
	t1 = [t45 * t43 * t62, t35, 0, 0, 0, 0; (-t31 * t59 - (-t41 * t44 * t45 * t61 + (t43 - 0.1e1) * t48 * t40) * t48 * t66) * t30, (t50 * t31 - (-t53 * t64 + t41 * t48 + (-t41 * t59 + t64) * t35) * t48 * t32) * t54 * t30, 0, 0, 0, 0; ((-t50 * t60 - t56) * t36 - (-t50 * t58 + t57) * t65) * t34, (-t36 * t47 + t49 * t65) * t34 * t62, 0, 0, t55 * t34, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:26:57
	% EndTime: 2019-10-10 09:26:57
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (640->28), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->42)
	t66 = qJ(2) + pkin(9);
	t62 = sin(t66);
	t82 = t62 ^ 2;
	t64 = cos(t66);
	t65 = pkin(10) + qJ(5);
	t63 = cos(t65);
	t69 = cos(qJ(1));
	t71 = t69 * t63;
	t61 = sin(t65);
	t68 = sin(qJ(1));
	t74 = t68 * t61;
	t49 = t64 * t74 + t71;
	t76 = t62 * t61;
	t45 = atan2(-t49, t76);
	t42 = sin(t45);
	t43 = cos(t45);
	t41 = -t42 * t49 + t43 * t76;
	t40 = 0.1e1 / t41 ^ 2;
	t72 = t69 * t61;
	t73 = t68 * t63;
	t52 = t64 * t72 - t73;
	t81 = t40 * t52;
	t79 = t43 * t49;
	t78 = t52 ^ 2 * t40;
	t56 = 0.1e1 / t61;
	t59 = 0.1e1 / t62;
	t77 = t56 * t59;
	t75 = t62 * t69;
	t53 = t64 * t71 + t74;
	t48 = 0.1e1 / t53 ^ 2;
	t70 = t69 ^ 2 * t82 * t48;
	t60 = 0.1e1 / t82;
	t57 = 0.1e1 / t61 ^ 2;
	t51 = t64 * t73 - t72;
	t47 = 0.1e1 / t53;
	t46 = 0.1e1 / (0.1e1 + t70);
	t44 = 0.1e1 / (t49 ^ 2 * t60 * t57 + 0.1e1);
	t39 = 0.1e1 / t41;
	t38 = (t49 * t56 * t60 * t64 + t68) * t44;
	t37 = 0.1e1 / (0.1e1 + t78);
	t36 = (t49 * t57 * t63 - t51 * t56) * t59 * t44;
	t1 = [-t52 * t44 * t77, t38, 0, 0, t36, 0; (-t49 * t39 - (-t42 + (t77 * t79 + t42) * t44) * t78) * t37, (t38 * t79 * t81 + (-t39 * t75 - (t43 * t64 + (-t38 + t68) * t62 * t42) * t81) * t61) * t37, 0, 0, (t53 * t39 - (t43 * t62 * t63 - t42 * t51 + (-t42 * t76 - t79) * t36) * t81) * t37, 0; (-t48 * t51 * t69 + t47 * t68) * t62 * t46, (-t47 * t64 * t69 - t63 * t70) * t46, 0, 0, -t52 * t48 * t46 * t75, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end