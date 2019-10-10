% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP4
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
%   Wie in S6RPRRPP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:14
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:55
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (238->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
	t43 = pkin(9) + qJ(3);
	t42 = cos(t43);
	t41 = sin(t43);
	t45 = sin(qJ(1));
	t53 = t45 * t41;
	t36 = atan2(-t53, -t42);
	t32 = sin(t36);
	t33 = cos(t36);
	t27 = -t32 * t53 - t33 * t42;
	t26 = 0.1e1 / t27 ^ 2;
	t47 = cos(qJ(1));
	t59 = t26 * t47 ^ 2;
	t46 = cos(qJ(4));
	t49 = t47 * t46;
	t44 = sin(qJ(4));
	t52 = t45 * t44;
	t35 = t42 * t49 + t52;
	t31 = 0.1e1 / t35 ^ 2;
	t50 = t47 * t44;
	t51 = t45 * t46;
	t34 = t42 * t50 - t51;
	t58 = t31 * t34;
	t57 = t32 * t42;
	t38 = t41 ^ 2;
	t56 = t38 / t42 ^ 2;
	t55 = t41 * t47;
	t37 = 0.1e1 / (t45 ^ 2 * t56 + 0.1e1);
	t54 = t45 * t37;
	t48 = t34 ^ 2 * t31 + 0.1e1;
	t39 = 0.1e1 / t42;
	t30 = 0.1e1 / t35;
	t29 = 0.1e1 / t48;
	t28 = (0.1e1 + t56) * t54;
	t25 = 0.1e1 / t27;
	t24 = 0.1e1 / (t38 * t59 + 0.1e1);
	t1 = [t39 * t37 * t55, 0, t28, 0, 0, 0; (-t25 * t53 - (-t33 * t38 * t39 * t54 + (t37 - 0.1e1) * t41 * t32) * t41 * t59) * t24, 0, (t42 * t25 - (-t45 * t57 + t33 * t41 + (-t33 * t53 + t57) * t28) * t41 * t26) * t47 * t24, 0, 0, 0; ((-t42 * t52 - t49) * t30 - (-t42 * t51 + t50) * t58) * t29, 0, (-t30 * t44 + t46 * t58) * t29 * t55, t48 * t29, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:55
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (286->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t52 = pkin(9) + qJ(3);
	t50 = cos(t52);
	t48 = sin(t52);
	t54 = sin(qJ(1));
	t61 = t54 * t48;
	t43 = atan2(-t61, -t50);
	t41 = sin(t43);
	t42 = cos(t43);
	t34 = -t41 * t61 - t42 * t50;
	t33 = 0.1e1 / t34 ^ 2;
	t55 = cos(qJ(1));
	t67 = t33 * t55 ^ 2;
	t53 = qJ(4) + pkin(10);
	t51 = cos(t53);
	t57 = t55 * t51;
	t49 = sin(t53);
	t60 = t54 * t49;
	t40 = t50 * t57 + t60;
	t38 = 0.1e1 / t40 ^ 2;
	t58 = t55 * t49;
	t59 = t54 * t51;
	t39 = t50 * t58 - t59;
	t66 = t38 * t39;
	t65 = t41 * t50;
	t45 = t48 ^ 2;
	t64 = t45 / t50 ^ 2;
	t63 = t48 * t55;
	t44 = 0.1e1 / (t54 ^ 2 * t64 + 0.1e1);
	t62 = t54 * t44;
	t56 = t39 ^ 2 * t38 + 0.1e1;
	t46 = 0.1e1 / t50;
	t37 = 0.1e1 / t40;
	t36 = (0.1e1 + t64) * t62;
	t35 = 0.1e1 / t56;
	t32 = 0.1e1 / t34;
	t31 = 0.1e1 / (t45 * t67 + 0.1e1);
	t1 = [t46 * t44 * t63, 0, t36, 0, 0, 0; (-t32 * t61 - (-t42 * t45 * t46 * t62 + (t44 - 0.1e1) * t48 * t41) * t48 * t67) * t31, 0, (t50 * t32 - (-t54 * t65 + t42 * t48 + (-t42 * t61 + t65) * t36) * t48 * t33) * t55 * t31, 0, 0, 0; ((-t50 * t60 - t57) * t37 - (-t50 * t59 + t58) * t66) * t35, 0, (-t37 * t49 + t51 * t66) * t35 * t63, t56 * t35, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:14:55
	% EndTime: 2019-10-10 01:14:56
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (640->28), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->42)
	t65 = pkin(9) + qJ(3);
	t61 = sin(t65);
	t82 = t61 ^ 2;
	t63 = cos(t65);
	t66 = qJ(4) + pkin(10);
	t64 = cos(t66);
	t69 = cos(qJ(1));
	t71 = t69 * t64;
	t62 = sin(t66);
	t68 = sin(qJ(1));
	t74 = t68 * t62;
	t49 = t63 * t74 + t71;
	t76 = t61 * t62;
	t45 = atan2(-t49, t76);
	t42 = sin(t45);
	t43 = cos(t45);
	t41 = -t42 * t49 + t43 * t76;
	t40 = 0.1e1 / t41 ^ 2;
	t72 = t69 * t62;
	t73 = t68 * t64;
	t52 = t63 * t72 - t73;
	t81 = t40 * t52;
	t79 = t43 * t49;
	t78 = t52 ^ 2 * t40;
	t57 = 0.1e1 / t61;
	t59 = 0.1e1 / t62;
	t77 = t57 * t59;
	t75 = t61 * t69;
	t53 = t63 * t71 + t74;
	t48 = 0.1e1 / t53 ^ 2;
	t70 = t69 ^ 2 * t82 * t48;
	t60 = 0.1e1 / t62 ^ 2;
	t58 = 0.1e1 / t82;
	t51 = t63 * t73 - t72;
	t47 = 0.1e1 / t53;
	t46 = 0.1e1 / (0.1e1 + t70);
	t44 = 0.1e1 / (t49 ^ 2 * t58 * t60 + 0.1e1);
	t39 = 0.1e1 / t41;
	t38 = (t49 * t58 * t59 * t63 + t68) * t44;
	t37 = 0.1e1 / (0.1e1 + t78);
	t36 = (t49 * t60 * t64 - t51 * t59) * t57 * t44;
	t1 = [-t52 * t44 * t77, 0, t38, t36, 0, 0; (-t49 * t39 - (-t42 + (t77 * t79 + t42) * t44) * t78) * t37, 0, (t38 * t79 * t81 + (-t39 * t75 - (t43 * t63 + (-t38 + t68) * t61 * t42) * t81) * t62) * t37, (t53 * t39 - (t43 * t61 * t64 - t42 * t51 + (-t42 * t76 - t79) * t36) * t81) * t37, 0, 0; (-t48 * t51 * t69 + t47 * t68) * t61 * t46, 0, (-t47 * t63 * t69 - t64 * t70) * t46, -t52 * t48 * t46 * t75, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end