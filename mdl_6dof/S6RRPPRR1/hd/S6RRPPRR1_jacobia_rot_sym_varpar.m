% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR1
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
%   Wie in S6RRPPRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (193->17), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->27)
	t36 = cos(qJ(1));
	t37 = t36 ^ 2;
	t32 = qJ(2) + pkin(10);
	t31 = cos(t32);
	t30 = sin(t32);
	t35 = sin(qJ(1));
	t38 = t35 * t30;
	t24 = atan2(-t38, -t31);
	t22 = sin(t24);
	t23 = cos(t24);
	t20 = -t22 * t38 - t23 * t31;
	t19 = 0.1e1 / t20 ^ 2;
	t43 = t19 * t30;
	t42 = t22 * t31;
	t27 = t30 ^ 2;
	t29 = 0.1e1 / t31 ^ 2;
	t41 = t27 * t29;
	t33 = t35 ^ 2;
	t40 = t33 / t37;
	t25 = 0.1e1 / (t33 * t41 + 0.1e1);
	t39 = t35 * t25;
	t28 = 0.1e1 / t31;
	t26 = 0.1e1 / (t29 * t40 + 0.1e1);
	t21 = (0.1e1 + t41) * t39;
	t18 = 0.1e1 / t20;
	t17 = 0.1e1 / (t37 * t27 * t19 + 0.1e1);
	t1 = [t36 * t30 * t28 * t25, t21, 0, 0, 0, 0; (-t18 * t38 - (-t23 * t27 * t28 * t39 + (t25 - 0.1e1) * t30 * t22) * t37 * t43) * t17, (t31 * t18 - (-t35 * t42 + t23 * t30 + (-t23 * t38 + t42) * t21) * t43) * t36 * t17, 0, 0, 0, 0; (-0.1e1 - t40) * t28 * t26, -0.1e1 / t36 * t29 * t26 * t38, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (737->24), mult. (876->55), div. (85->9), fcn. (1282->11), ass. (0->38)
	t87 = qJ(2) + pkin(10);
	t73 = cos(t87);
	t75 = sin(qJ(5));
	t84 = sin(t87);
	t95 = cos(qJ(5));
	t67 = -t73 * t75 + t84 * t95;
	t76 = sin(qJ(1));
	t59 = t67 * t76;
	t66 = t73 * t95 + t84 * t75;
	t54 = atan2(t59, t66);
	t51 = sin(t54);
	t52 = cos(t54);
	t49 = t51 * t59 + t52 * t66;
	t48 = 0.1e1 / t49 ^ 2;
	t96 = cos(qJ(1));
	t63 = t67 * t96;
	t89 = t63 ^ 2 * t48;
	t46 = 0.1e1 / (0.1e1 + t89);
	t47 = 0.1e1 / t49;
	t61 = t66 * t76;
	t62 = t66 * t96;
	t92 = t52 * t59;
	t65 = 0.1e1 / t66 ^ 2;
	t53 = 0.1e1 / (t59 ^ 2 * t65 + 0.1e1);
	t64 = 0.1e1 / t66;
	t98 = (-t59 * t65 * t67 - t61 * t64) * t53;
	t101 = (t48 * t63 * (-t51 * t61 + t52 * t67 + (-t51 * t66 + t92) * t98) + t62 * t47) * t46;
	t74 = sin(qJ(6));
	t77 = cos(qJ(6));
	t58 = t62 * t77 - t76 * t74;
	t56 = 0.1e1 / t58 ^ 2;
	t57 = t62 * t74 + t76 * t77;
	t86 = t57 ^ 2 * t56 + 0.1e1;
	t50 = 0.1e1 / t86;
	t55 = 0.1e1 / t58;
	t91 = t56 * t57;
	t97 = (-t74 * t55 + t77 * t91) * t50 * t63;
	t1 = [t63 * t64 * t53, -t98, 0, 0, t98, 0; (t59 * t47 - (-t51 + (-t64 * t92 + t51) * t53) * t89) * t46, -t101, 0, 0, t101, 0; ((-t61 * t74 + t96 * t77) * t55 - (-t61 * t77 - t96 * t74) * t91) * t50, t97, 0, 0, -t97, t86 * t50;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end