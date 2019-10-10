% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR7
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
%   Wie in S6RPPRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (198->20), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->35)
	t35 = pkin(9) + qJ(4);
	t33 = sin(t35);
	t34 = cos(t35);
	t39 = cos(qJ(1));
	t42 = t39 * t34;
	t28 = atan2(-t42, t33);
	t26 = sin(t28);
	t27 = cos(t28);
	t19 = -t26 * t42 + t27 * t33;
	t18 = 0.1e1 / t19 ^ 2;
	t38 = sin(qJ(1));
	t50 = t18 * t38 ^ 2;
	t36 = sin(pkin(10));
	t41 = t39 * t36;
	t37 = cos(pkin(10));
	t44 = t38 * t37;
	t25 = t33 * t44 + t41;
	t23 = 0.1e1 / t25 ^ 2;
	t40 = t39 * t37;
	t45 = t38 * t36;
	t24 = t33 * t45 - t40;
	t49 = t23 * t24;
	t48 = t26 * t33;
	t32 = t34 ^ 2;
	t47 = 0.1e1 / t33 ^ 2 * t32;
	t46 = t34 * t38;
	t29 = 0.1e1 / (t39 ^ 2 * t47 + 0.1e1);
	t43 = t39 * t29;
	t30 = 0.1e1 / t33;
	t22 = 0.1e1 / t25;
	t21 = 0.1e1 / (t24 ^ 2 * t23 + 0.1e1);
	t20 = (0.1e1 + t47) * t43;
	t17 = 0.1e1 / t19;
	t16 = 0.1e1 / (t32 * t50 + 0.1e1);
	t1 = [t30 * t29 * t46, 0, 0, t20, 0, 0; (-t17 * t42 + (-t27 * t30 * t32 * t43 + (-t29 + 0.1e1) * t34 * t26) * t34 * t50) * t16, 0, 0, (t33 * t17 + (t39 * t48 + t27 * t34 + (-t27 * t42 - t48) * t20) * t34 * t18) * t38 * t16, 0, 0; ((t33 * t41 + t44) * t22 - (t33 * t40 - t45) * t49) * t21, 0, 0, (t22 * t36 - t37 * t49) * t21 * t46, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:44:24
	% EndTime: 2019-10-09 23:44:24
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (263->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t50 = pkin(9) + qJ(4);
	t46 = sin(t50);
	t48 = cos(t50);
	t52 = cos(qJ(1));
	t54 = t52 * t48;
	t40 = atan2(-t54, t46);
	t38 = sin(t40);
	t39 = cos(t40);
	t31 = -t38 * t54 + t39 * t46;
	t30 = 0.1e1 / t31 ^ 2;
	t51 = sin(qJ(1));
	t64 = t30 * t51 ^ 2;
	t49 = pkin(10) + qJ(6);
	t45 = sin(t49);
	t56 = t52 * t45;
	t47 = cos(t49);
	t58 = t51 * t47;
	t37 = t46 * t58 + t56;
	t35 = 0.1e1 / t37 ^ 2;
	t55 = t52 * t47;
	t59 = t51 * t45;
	t36 = t46 * t59 - t55;
	t63 = t35 * t36;
	t62 = t38 * t46;
	t44 = t48 ^ 2;
	t61 = 0.1e1 / t46 ^ 2 * t44;
	t60 = t48 * t51;
	t41 = 0.1e1 / (t52 ^ 2 * t61 + 0.1e1);
	t57 = t52 * t41;
	t53 = t36 ^ 2 * t35 + 0.1e1;
	t42 = 0.1e1 / t46;
	t34 = 0.1e1 / t37;
	t33 = (0.1e1 + t61) * t57;
	t32 = 0.1e1 / t53;
	t29 = 0.1e1 / t31;
	t28 = 0.1e1 / (t44 * t64 + 0.1e1);
	t1 = [t42 * t41 * t60, 0, 0, t33, 0, 0; (-t29 * t54 + (-t39 * t42 * t44 * t57 + (-t41 + 0.1e1) * t48 * t38) * t48 * t64) * t28, 0, 0, (t46 * t29 + (t52 * t62 + t39 * t48 + (-t39 * t54 - t62) * t33) * t48 * t30) * t51 * t28, 0, 0; ((t46 * t56 + t58) * t34 - (t46 * t55 - t59) * t63) * t32, 0, 0, (t34 * t45 - t47 * t63) * t32 * t60, 0, t53 * t32;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end