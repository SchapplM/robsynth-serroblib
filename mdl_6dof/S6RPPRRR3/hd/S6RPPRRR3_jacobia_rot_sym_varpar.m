% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR3
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
%   Wie in S6RPPRRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRR3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (195->20), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
	t36 = qJ(1) + pkin(10);
	t34 = sin(t36);
	t54 = t34 ^ 2;
	t41 = sin(qJ(4));
	t35 = cos(t36);
	t43 = cos(qJ(4));
	t48 = t35 * t43;
	t32 = atan2(-t48, t41);
	t30 = sin(t32);
	t31 = cos(t32);
	t24 = -t30 * t48 + t31 * t41;
	t22 = 0.1e1 / t24 ^ 2;
	t53 = t22 * t43;
	t40 = sin(qJ(5));
	t42 = cos(qJ(5));
	t45 = t41 * t42;
	t29 = t34 * t45 + t35 * t40;
	t27 = 0.1e1 / t29 ^ 2;
	t46 = t40 * t41;
	t28 = t34 * t46 - t35 * t42;
	t52 = t27 * t28;
	t51 = t30 * t41;
	t50 = t34 * t43;
	t39 = t43 ^ 2;
	t47 = 0.1e1 / t41 ^ 2 * t39;
	t33 = 0.1e1 / (t35 ^ 2 * t47 + 0.1e1);
	t49 = t35 * t33;
	t44 = t28 ^ 2 * t27 + 0.1e1;
	t37 = 0.1e1 / t41;
	t26 = 0.1e1 / t29;
	t25 = (0.1e1 + t47) * t49;
	t23 = 0.1e1 / t44;
	t21 = 0.1e1 / t24;
	t20 = 0.1e1 / (t54 * t39 * t22 + 0.1e1);
	t1 = [t37 * t33 * t50, 0, 0, t25, 0, 0; (-t21 * t48 + (-t31 * t37 * t39 * t49 + (-t33 + 0.1e1) * t43 * t30) * t54 * t53) * t20, 0, 0, (t41 * t21 + (t35 * t51 + t31 * t43 + (-t31 * t48 - t51) * t25) * t53) * t34 * t20, 0, 0; ((t34 * t42 + t35 * t46) * t26 - (-t34 * t40 + t35 * t45) * t52) * t23, 0, 0, (t26 * t40 - t42 * t52) * t23 * t50, t44 * t23, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:04:49
	% EndTime: 2019-10-10 00:04:49
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (281->21), mult. (251->57), div. (57->9), fcn. (367->9), ass. (0->37)
	t53 = qJ(1) + pkin(10);
	t49 = sin(t53);
	t70 = t49 ^ 2;
	t58 = sin(qJ(4));
	t50 = cos(t53);
	t59 = cos(qJ(4));
	t64 = t50 * t59;
	t47 = atan2(-t64, t58);
	t45 = sin(t47);
	t46 = cos(t47);
	t39 = -t45 * t64 + t46 * t58;
	t38 = 0.1e1 / t39 ^ 2;
	t69 = t38 * t59;
	t57 = qJ(5) + qJ(6);
	t51 = sin(t57);
	t52 = cos(t57);
	t62 = t52 * t58;
	t44 = t49 * t62 + t50 * t51;
	t42 = 0.1e1 / t44 ^ 2;
	t63 = t51 * t58;
	t43 = t49 * t63 - t50 * t52;
	t68 = t42 * t43;
	t67 = t45 * t58;
	t66 = t49 * t59;
	t56 = t59 ^ 2;
	t61 = 0.1e1 / t58 ^ 2 * t56;
	t48 = 0.1e1 / (t50 ^ 2 * t61 + 0.1e1);
	t65 = t50 * t48;
	t60 = t43 ^ 2 * t42 + 0.1e1;
	t54 = 0.1e1 / t58;
	t41 = 0.1e1 / t44;
	t40 = (0.1e1 + t61) * t65;
	t37 = 0.1e1 / t39;
	t36 = 0.1e1 / t60;
	t35 = 0.1e1 / (t70 * t56 * t38 + 0.1e1);
	t34 = t60 * t36;
	t1 = [t54 * t48 * t66, 0, 0, t40, 0, 0; (-t37 * t64 + (-t46 * t54 * t56 * t65 + (-t48 + 0.1e1) * t59 * t45) * t70 * t69) * t35, 0, 0, (t58 * t37 + (t50 * t67 + t46 * t59 + (-t46 * t64 - t67) * t40) * t69) * t49 * t35, 0, 0; ((t49 * t52 + t50 * t63) * t41 - (-t49 * t51 + t50 * t62) * t68) * t36, 0, 0, (t41 * t51 - t52 * t68) * t36 * t66, t34, t34;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end