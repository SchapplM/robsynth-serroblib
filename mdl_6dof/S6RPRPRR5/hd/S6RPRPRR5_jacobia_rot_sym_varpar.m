% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR5
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
%   Wie in S6RPRPRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:53
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (193->17), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->27)
	t33 = cos(qJ(1));
	t34 = t33 ^ 2;
	t29 = pkin(10) + qJ(3);
	t28 = cos(t29);
	t27 = sin(t29);
	t32 = sin(qJ(1));
	t35 = t32 * t27;
	t21 = atan2(-t35, -t28);
	t19 = sin(t21);
	t20 = cos(t21);
	t17 = -t19 * t35 - t20 * t28;
	t16 = 0.1e1 / t17 ^ 2;
	t40 = t16 * t27;
	t39 = t19 * t28;
	t24 = t27 ^ 2;
	t26 = 0.1e1 / t28 ^ 2;
	t38 = t24 * t26;
	t30 = t32 ^ 2;
	t37 = t30 / t34;
	t22 = 0.1e1 / (t30 * t38 + 0.1e1);
	t36 = t32 * t22;
	t25 = 0.1e1 / t28;
	t23 = 0.1e1 / (t26 * t37 + 0.1e1);
	t18 = (0.1e1 + t38) * t36;
	t15 = 0.1e1 / t17;
	t14 = 0.1e1 / (t34 * t24 * t16 + 0.1e1);
	t1 = [t33 * t27 * t25 * t22, 0, t18, 0, 0, 0; (-t15 * t35 - (-t20 * t24 * t25 * t36 + (t22 - 0.1e1) * t27 * t19) * t34 * t40) * t14, 0, (t28 * t15 - (-t32 * t39 + t20 * t27 + (-t20 * t35 + t39) * t18) * t40) * t33 * t14, 0, 0, 0; (-0.1e1 - t37) * t25 * t23, 0, -0.1e1 / t33 * t26 * t23 * t35, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:53:10
	% EndTime: 2019-10-10 00:53:10
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (737->24), mult. (876->55), div. (85->9), fcn. (1282->11), ass. (0->38)
	t84 = pkin(10) + qJ(3);
	t70 = cos(t84);
	t72 = sin(qJ(5));
	t81 = sin(t84);
	t92 = cos(qJ(5));
	t64 = -t70 * t72 + t81 * t92;
	t73 = sin(qJ(1));
	t56 = t64 * t73;
	t63 = t70 * t92 + t72 * t81;
	t51 = atan2(t56, t63);
	t48 = sin(t51);
	t49 = cos(t51);
	t46 = t48 * t56 + t49 * t63;
	t45 = 0.1e1 / t46 ^ 2;
	t93 = cos(qJ(1));
	t60 = t64 * t93;
	t86 = t60 ^ 2 * t45;
	t43 = 0.1e1 / (0.1e1 + t86);
	t44 = 0.1e1 / t46;
	t58 = t63 * t73;
	t59 = t63 * t93;
	t89 = t49 * t56;
	t62 = 0.1e1 / t63 ^ 2;
	t50 = 0.1e1 / (t56 ^ 2 * t62 + 0.1e1);
	t61 = 0.1e1 / t63;
	t95 = (-t56 * t62 * t64 - t58 * t61) * t50;
	t98 = (t45 * t60 * (-t48 * t58 + t49 * t64 + (-t48 * t63 + t89) * t95) + t59 * t44) * t43;
	t71 = sin(qJ(6));
	t74 = cos(qJ(6));
	t55 = t59 * t74 - t71 * t73;
	t53 = 0.1e1 / t55 ^ 2;
	t54 = t59 * t71 + t73 * t74;
	t83 = t53 * t54 ^ 2 + 0.1e1;
	t47 = 0.1e1 / t83;
	t52 = 0.1e1 / t55;
	t88 = t53 * t54;
	t94 = (-t52 * t71 + t74 * t88) * t47 * t60;
	t1 = [t60 * t61 * t50, 0, -t95, 0, t95, 0; (t56 * t44 - (-t48 + (-t61 * t89 + t48) * t50) * t86) * t43, 0, -t98, 0, t98, 0; ((-t58 * t71 + t74 * t93) * t52 - (-t58 * t74 - t71 * t93) * t88) * t47, 0, t94, 0, -t94, t83 * t47;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end