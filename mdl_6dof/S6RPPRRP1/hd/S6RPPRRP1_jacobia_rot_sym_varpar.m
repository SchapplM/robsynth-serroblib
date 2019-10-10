% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRP1
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
%   Wie in S6RPPRRP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:47
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:47:46
	% EndTime: 2019-10-09 23:47:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:47:46
	% EndTime: 2019-10-09 23:47:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:47:46
	% EndTime: 2019-10-09 23:47:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:47:46
	% EndTime: 2019-10-09 23:47:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:47:46
	% EndTime: 2019-10-09 23:47:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:47:46
	% EndTime: 2019-10-09 23:47:47
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (347->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t46 = pkin(10) + qJ(4);
	t44 = cos(t46);
	t42 = sin(t46);
	t47 = qJ(1) + pkin(9);
	t43 = sin(t47);
	t55 = t43 * t42;
	t37 = atan2(-t55, -t44);
	t35 = sin(t37);
	t36 = cos(t37);
	t28 = -t35 * t55 - t36 * t44;
	t27 = 0.1e1 / t28 ^ 2;
	t45 = cos(t47);
	t61 = t27 * t45 ^ 2;
	t49 = cos(qJ(5));
	t51 = t45 * t49;
	t48 = sin(qJ(5));
	t54 = t43 * t48;
	t34 = t44 * t51 + t54;
	t32 = 0.1e1 / t34 ^ 2;
	t52 = t45 * t48;
	t53 = t43 * t49;
	t33 = t44 * t52 - t53;
	t60 = t32 * t33;
	t59 = t35 * t44;
	t39 = t42 ^ 2;
	t58 = t39 / t44 ^ 2;
	t57 = t42 * t45;
	t38 = 0.1e1 / (t43 ^ 2 * t58 + 0.1e1);
	t56 = t43 * t38;
	t50 = t33 ^ 2 * t32 + 0.1e1;
	t40 = 0.1e1 / t44;
	t31 = 0.1e1 / t34;
	t30 = 0.1e1 / t50;
	t29 = (0.1e1 + t58) * t56;
	t26 = 0.1e1 / t28;
	t25 = 0.1e1 / (t39 * t61 + 0.1e1);
	t1 = [t40 * t38 * t57, 0, 0, t29, 0, 0; (-t26 * t55 - (-t36 * t39 * t40 * t56 + (t38 - 0.1e1) * t42 * t35) * t42 * t61) * t25, 0, 0, (t44 * t26 - (-t43 * t59 + t36 * t42 + (-t36 * t55 + t59) * t29) * t42 * t27) * t45 * t25, 0, 0; ((-t44 * t54 - t51) * t31 - (-t44 * t53 + t52) * t60) * t30, 0, 0, (-t31 * t48 + t49 * t60) * t30 * t57, t50 * t30, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:47:46
	% EndTime: 2019-10-09 23:47:47
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (347->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t50 = pkin(10) + qJ(4);
	t48 = cos(t50);
	t46 = sin(t50);
	t51 = qJ(1) + pkin(9);
	t47 = sin(t51);
	t59 = t47 * t46;
	t41 = atan2(-t59, -t48);
	t39 = sin(t41);
	t40 = cos(t41);
	t32 = -t39 * t59 - t40 * t48;
	t31 = 0.1e1 / t32 ^ 2;
	t49 = cos(t51);
	t65 = t31 * t49 ^ 2;
	t53 = cos(qJ(5));
	t55 = t49 * t53;
	t52 = sin(qJ(5));
	t58 = t47 * t52;
	t38 = t48 * t55 + t58;
	t36 = 0.1e1 / t38 ^ 2;
	t56 = t49 * t52;
	t57 = t47 * t53;
	t37 = t48 * t56 - t57;
	t64 = t36 * t37;
	t63 = t39 * t48;
	t43 = t46 ^ 2;
	t62 = t43 / t48 ^ 2;
	t61 = t46 * t49;
	t42 = 0.1e1 / (t47 ^ 2 * t62 + 0.1e1);
	t60 = t47 * t42;
	t54 = t37 ^ 2 * t36 + 0.1e1;
	t44 = 0.1e1 / t48;
	t35 = 0.1e1 / t38;
	t34 = 0.1e1 / t54;
	t33 = (0.1e1 + t62) * t60;
	t30 = 0.1e1 / t32;
	t29 = 0.1e1 / (t43 * t65 + 0.1e1);
	t1 = [t44 * t42 * t61, 0, 0, t33, 0, 0; (-t30 * t59 - (-t40 * t43 * t44 * t60 + (t42 - 0.1e1) * t46 * t39) * t46 * t65) * t29, 0, 0, (t48 * t30 - (-t47 * t63 + t40 * t46 + (-t40 * t59 + t63) * t33) * t46 * t31) * t49 * t29, 0, 0; ((-t48 * t58 - t55) * t35 - (-t48 * t57 + t56) * t64) * t34, 0, 0, (-t35 * t52 + t53 * t64) * t34 * t61, t54 * t34, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end