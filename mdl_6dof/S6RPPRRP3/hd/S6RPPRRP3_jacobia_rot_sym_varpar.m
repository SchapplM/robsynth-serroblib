% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRP3
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
%   Wie in S6RPPRRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:51
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (195->20), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
	t36 = qJ(1) + pkin(9);
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
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (379->25), mult. (509->71), div. (100->11), fcn. (770->9), ass. (0->39)
	t55 = cos(qJ(4));
	t66 = t55 ^ 2;
	t46 = qJ(1) + pkin(9);
	t44 = sin(t46);
	t45 = cos(t46);
	t54 = cos(qJ(5));
	t52 = sin(qJ(5));
	t53 = sin(qJ(4));
	t59 = t52 * t53;
	t40 = t44 * t54 + t45 * t59;
	t57 = t55 * t52;
	t35 = atan2(t40, t57);
	t32 = sin(t35);
	t33 = cos(t35);
	t30 = t32 * t40 + t33 * t57;
	t29 = 0.1e1 / t30 ^ 2;
	t38 = t44 * t59 - t45 * t54;
	t65 = t29 * t38;
	t63 = t33 * t40;
	t62 = t38 ^ 2 * t29;
	t61 = t44 * t55;
	t47 = 0.1e1 / t52;
	t50 = 0.1e1 / t55;
	t60 = t47 * t50;
	t58 = t53 * t54;
	t39 = t44 * t58 + t45 * t52;
	t37 = 0.1e1 / t39 ^ 2;
	t56 = t44 ^ 2 * t66 * t37;
	t51 = 0.1e1 / t66;
	t48 = 0.1e1 / t52 ^ 2;
	t41 = -t44 * t52 + t45 * t58;
	t36 = 0.1e1 / t39;
	t34 = 0.1e1 / (t40 ^ 2 * t51 * t48 + 0.1e1);
	t31 = 0.1e1 / (0.1e1 + t56);
	t28 = 0.1e1 / t30;
	t27 = (t40 * t47 * t51 * t53 + t45) * t34;
	t26 = 0.1e1 / (0.1e1 + t62);
	t25 = (-t40 * t48 * t54 + t41 * t47) * t50 * t34;
	t1 = [-t38 * t34 * t60, 0, 0, t27, t25, 0; (t40 * t28 - (-t32 + (-t60 * t63 + t32) * t34) * t62) * t26, 0, 0, (-t27 * t63 * t65 + (t28 * t61 - (-t33 * t53 + (-t27 + t45) * t32 * t55) * t65) * t52) * t26, (t39 * t28 - (t33 * t55 * t54 + t32 * t41 + (-t32 * t57 + t63) * t25) * t65) * t26, 0; (-t37 * t41 * t44 + t36 * t45) * t55 * t31, 0, 0, (-t36 * t44 * t53 - t54 * t56) * t31, t38 * t37 * t31 * t61, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end