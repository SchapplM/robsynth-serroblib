% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR7
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
%   Wie in S6RRPRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:15
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (82->16), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->26)
	t29 = cos(qJ(1));
	t30 = t29 ^ 2;
	t28 = cos(qJ(2));
	t26 = sin(qJ(2));
	t27 = sin(qJ(1));
	t31 = t27 * t26;
	t18 = atan2(-t31, -t28);
	t16 = sin(t18);
	t17 = cos(t18);
	t14 = -t16 * t31 - t17 * t28;
	t13 = 0.1e1 / t14 ^ 2;
	t36 = t13 * t26;
	t35 = t16 * t28;
	t21 = t26 ^ 2;
	t24 = 0.1e1 / t28 ^ 2;
	t34 = t21 * t24;
	t22 = t27 ^ 2;
	t33 = t22 / t30;
	t19 = 0.1e1 / (t22 * t34 + 0.1e1);
	t32 = t27 * t19;
	t23 = 0.1e1 / t28;
	t20 = 0.1e1 / (t24 * t33 + 0.1e1);
	t15 = (0.1e1 + t34) * t32;
	t12 = 0.1e1 / t14;
	t11 = 0.1e1 / (t30 * t21 * t13 + 0.1e1);
	t1 = [t29 * t26 * t23 * t19, t15, 0, 0, 0, 0; (-t12 * t31 - (-t17 * t21 * t23 * t32 + (t19 - 0.1e1) * t26 * t16) * t30 * t36) * t11, (t28 * t12 - (-t27 * t35 + t17 * t26 + (-t17 * t31 + t35) * t15) * t36) * t29 * t11, 0, 0, 0, 0; (-0.1e1 - t33) * t23 * t20, -0.1e1 / t29 * t24 * t20 * t31, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (737->24), mult. (876->55), div. (85->9), fcn. (1282->11), ass. (0->38)
	t87 = qJ(4) + pkin(10);
	t73 = sin(t87);
	t77 = cos(qJ(2));
	t84 = cos(t87);
	t95 = sin(qJ(2));
	t67 = -t77 * t73 + t95 * t84;
	t75 = sin(qJ(1));
	t59 = t67 * t75;
	t66 = t95 * t73 + t77 * t84;
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
	t61 = t66 * t75;
	t62 = t66 * t96;
	t92 = t52 * t59;
	t65 = 0.1e1 / t66 ^ 2;
	t53 = 0.1e1 / (t59 ^ 2 * t65 + 0.1e1);
	t64 = 0.1e1 / t66;
	t98 = (-t59 * t65 * t67 - t61 * t64) * t53;
	t101 = (t48 * t63 * (-t51 * t61 + t52 * t67 + (-t51 * t66 + t92) * t98) + t62 * t47) * t46;
	t74 = sin(qJ(6));
	t76 = cos(qJ(6));
	t58 = t62 * t76 - t75 * t74;
	t56 = 0.1e1 / t58 ^ 2;
	t57 = t62 * t74 + t75 * t76;
	t86 = t57 ^ 2 * t56 + 0.1e1;
	t50 = 0.1e1 / t86;
	t55 = 0.1e1 / t58;
	t91 = t56 * t57;
	t97 = (-t74 * t55 + t76 * t91) * t50 * t63;
	t1 = [t63 * t64 * t53, -t98, 0, t98, 0, 0; (t59 * t47 - (-t51 + (-t64 * t92 + t51) * t53) * t89) * t46, -t101, 0, t101, 0, 0; ((-t61 * t74 + t96 * t76) * t55 - (-t61 * t76 - t96 * t74) * t91) * t50, t97, 0, -t97, 0, t86 * t50;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end