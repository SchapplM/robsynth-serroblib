% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP5
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
%   Wie in S6RPRPRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:36
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
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
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (286->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t49 = pkin(9) + qJ(3);
	t47 = cos(t49);
	t45 = sin(t49);
	t50 = sin(qJ(1));
	t56 = t50 * t45;
	t39 = atan2(-t56, -t47);
	t37 = sin(t39);
	t38 = cos(t39);
	t30 = -t37 * t56 - t38 * t47;
	t29 = 0.1e1 / t30 ^ 2;
	t51 = cos(qJ(1));
	t63 = t29 * t51 ^ 2;
	t48 = pkin(10) + qJ(5);
	t46 = cos(t48);
	t53 = t51 * t46;
	t44 = sin(t48);
	t57 = t50 * t44;
	t36 = t47 * t53 + t57;
	t34 = 0.1e1 / t36 ^ 2;
	t54 = t51 * t44;
	t55 = t50 * t46;
	t35 = t47 * t54 - t55;
	t62 = t34 * t35;
	t61 = t37 * t47;
	t41 = t45 ^ 2;
	t60 = t41 / t47 ^ 2;
	t59 = t45 * t51;
	t40 = 0.1e1 / (t50 ^ 2 * t60 + 0.1e1);
	t58 = t50 * t40;
	t52 = t35 ^ 2 * t34 + 0.1e1;
	t42 = 0.1e1 / t47;
	t33 = 0.1e1 / t36;
	t32 = (0.1e1 + t60) * t58;
	t31 = 0.1e1 / t52;
	t28 = 0.1e1 / t30;
	t27 = 0.1e1 / (t41 * t63 + 0.1e1);
	t1 = [t42 * t40 * t59, 0, t32, 0, 0, 0; (-t28 * t56 - (-t38 * t41 * t42 * t58 + (t40 - 0.1e1) * t45 * t37) * t45 * t63) * t27, 0, (t47 * t28 - (-t50 * t61 + t38 * t45 + (-t38 * t56 + t61) * t32) * t45 * t29) * t51 * t27, 0, 0, 0; ((-t47 * t57 - t53) * t33 - (-t47 * t55 + t54) * t62) * t31, 0, (-t33 * t44 + t46 * t62) * t31 * t59, 0, t52 * t31, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:35:58
	% EndTime: 2019-10-10 00:35:58
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (640->28), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->42)
	t62 = pkin(9) + qJ(3);
	t58 = sin(t62);
	t78 = t58 ^ 2;
	t60 = cos(t62);
	t61 = pkin(10) + qJ(5);
	t59 = cos(t61);
	t65 = cos(qJ(1));
	t67 = t65 * t59;
	t57 = sin(t61);
	t64 = sin(qJ(1));
	t70 = t64 * t57;
	t45 = t60 * t70 + t67;
	t72 = t58 * t57;
	t41 = atan2(-t45, t72);
	t38 = sin(t41);
	t39 = cos(t41);
	t37 = -t38 * t45 + t39 * t72;
	t36 = 0.1e1 / t37 ^ 2;
	t68 = t65 * t57;
	t69 = t64 * t59;
	t48 = t60 * t68 - t69;
	t77 = t36 * t48;
	t75 = t39 * t45;
	t74 = t48 ^ 2 * t36;
	t52 = 0.1e1 / t57;
	t55 = 0.1e1 / t58;
	t73 = t52 * t55;
	t71 = t58 * t65;
	t49 = t60 * t67 + t70;
	t44 = 0.1e1 / t49 ^ 2;
	t66 = t65 ^ 2 * t78 * t44;
	t56 = 0.1e1 / t78;
	t53 = 0.1e1 / t57 ^ 2;
	t47 = t60 * t69 - t68;
	t43 = 0.1e1 / t49;
	t42 = 0.1e1 / (0.1e1 + t66);
	t40 = 0.1e1 / (t45 ^ 2 * t56 * t53 + 0.1e1);
	t35 = 0.1e1 / t37;
	t34 = (t45 * t52 * t56 * t60 + t64) * t40;
	t33 = 0.1e1 / (0.1e1 + t74);
	t32 = (t45 * t53 * t59 - t47 * t52) * t55 * t40;
	t1 = [-t48 * t40 * t73, 0, t34, 0, t32, 0; (-t45 * t35 - (-t38 + (t73 * t75 + t38) * t40) * t74) * t33, 0, (t34 * t75 * t77 + (-t35 * t71 - (t39 * t60 + (-t34 + t64) * t58 * t38) * t77) * t57) * t33, 0, (t49 * t35 - (t39 * t58 * t59 - t38 * t47 + (-t38 * t72 - t75) * t32) * t77) * t33, 0; (-t44 * t47 * t65 + t43 * t64) * t58 * t42, 0, (-t43 * t60 * t65 - t59 * t66) * t42, 0, -t48 * t44 * t42 * t71, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end