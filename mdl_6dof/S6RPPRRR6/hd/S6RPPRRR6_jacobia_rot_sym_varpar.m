% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR6
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
%   Wie in S6RPPRRR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (63->18), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->35)
	t40 = sin(qJ(4));
	t41 = sin(qJ(1));
	t43 = cos(qJ(4));
	t49 = t41 * t43;
	t35 = atan2(t49, t40);
	t32 = sin(t35);
	t33 = cos(t35);
	t25 = t32 * t49 + t33 * t40;
	t24 = 0.1e1 / t25 ^ 2;
	t44 = cos(qJ(1));
	t56 = t24 * t44 ^ 2;
	t42 = cos(qJ(5));
	t46 = t44 * t42;
	t39 = sin(qJ(5));
	t51 = t41 * t39;
	t31 = t40 * t46 - t51;
	t29 = 0.1e1 / t31 ^ 2;
	t47 = t44 * t39;
	t50 = t41 * t42;
	t30 = t40 * t47 + t50;
	t55 = t29 * t30;
	t54 = t32 * t40;
	t38 = t43 ^ 2;
	t53 = 0.1e1 / t40 ^ 2 * t38;
	t34 = 0.1e1 / (t41 ^ 2 * t53 + 0.1e1);
	t52 = t41 * t34;
	t48 = t43 * t44;
	t45 = t30 ^ 2 * t29 + 0.1e1;
	t36 = 0.1e1 / t40;
	t28 = 0.1e1 / t31;
	t27 = (-0.1e1 - t53) * t52;
	t26 = 0.1e1 / t45;
	t23 = 0.1e1 / t25;
	t22 = 0.1e1 / (t38 * t56 + 0.1e1);
	t1 = [t36 * t34 * t48, 0, 0, t27, 0, 0; (t23 * t49 + (t33 * t36 * t38 * t52 + (-t34 + 0.1e1) * t43 * t32) * t43 * t56) * t22, 0, 0, (t40 * t23 + (-t41 * t54 + t33 * t43 + (t33 * t49 - t54) * t27) * t43 * t24) * t44 * t22, 0, 0; ((-t40 * t51 + t46) * t28 - (-t40 * t50 - t47) * t55) * t26, 0, 0, (t28 * t39 - t42 * t55) * t26 * t48, t45 * t26, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:10:02
	% EndTime: 2019-10-10 00:10:02
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (135->19), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->37)
	t54 = sin(qJ(4));
	t55 = sin(qJ(1));
	t56 = cos(qJ(4));
	t62 = t55 * t56;
	t47 = atan2(t62, t54);
	t44 = sin(t47);
	t45 = cos(t47);
	t38 = t44 * t62 + t45 * t54;
	t37 = 0.1e1 / t38 ^ 2;
	t57 = cos(qJ(1));
	t69 = t37 * t57 ^ 2;
	t53 = qJ(5) + qJ(6);
	t49 = cos(t53);
	t59 = t57 * t49;
	t48 = sin(t53);
	t64 = t55 * t48;
	t43 = t54 * t59 - t64;
	t41 = 0.1e1 / t43 ^ 2;
	t60 = t57 * t48;
	t63 = t55 * t49;
	t42 = t54 * t60 + t63;
	t68 = t41 * t42;
	t67 = t44 * t54;
	t52 = t56 ^ 2;
	t66 = 0.1e1 / t54 ^ 2 * t52;
	t46 = 0.1e1 / (t55 ^ 2 * t66 + 0.1e1);
	t65 = t55 * t46;
	t61 = t56 * t57;
	t58 = t42 ^ 2 * t41 + 0.1e1;
	t50 = 0.1e1 / t54;
	t40 = 0.1e1 / t43;
	t39 = (-0.1e1 - t66) * t65;
	t36 = 0.1e1 / t38;
	t35 = 0.1e1 / t58;
	t34 = 0.1e1 / (t52 * t69 + 0.1e1);
	t33 = t58 * t35;
	t1 = [t50 * t46 * t61, 0, 0, t39, 0, 0; (t36 * t62 + (t45 * t50 * t52 * t65 + (-t46 + 0.1e1) * t56 * t44) * t56 * t69) * t34, 0, 0, (t54 * t36 + (-t55 * t67 + t45 * t56 + (t45 * t62 - t67) * t39) * t56 * t37) * t57 * t34, 0, 0; ((-t54 * t64 + t59) * t40 - (-t54 * t63 - t60) * t68) * t35, 0, 0, (t40 * t48 - t49 * t68) * t35 * t61, t33, t33;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end