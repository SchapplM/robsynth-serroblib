% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR2
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
%   Wie in S6RPRPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:17
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (263->18), mult. (149->39), div. (43->9), fcn. (227->7), ass. (0->29)
	t38 = qJ(1) + pkin(9);
	t36 = cos(t38);
	t32 = t36 ^ 2;
	t37 = qJ(3) + pkin(10);
	t35 = cos(t37);
	t33 = sin(t37);
	t34 = sin(t38);
	t41 = t34 * t33;
	t25 = atan2(-t41, -t35);
	t23 = sin(t25);
	t24 = cos(t25);
	t21 = -t23 * t41 - t24 * t35;
	t20 = 0.1e1 / t21 ^ 2;
	t47 = t20 * t33;
	t46 = t23 * t35;
	t28 = t33 ^ 2;
	t40 = t35 ^ 2;
	t45 = t28 / t40;
	t39 = t34 ^ 2;
	t44 = 0.1e1 / t39 * t32;
	t43 = t33 * t36;
	t26 = 0.1e1 / (t39 * t45 + 0.1e1);
	t42 = t34 * t26;
	t30 = 0.1e1 / t35;
	t27 = 0.1e1 / (t40 * t44 + 0.1e1);
	t22 = (0.1e1 + t45) * t42;
	t19 = 0.1e1 / t21;
	t18 = 0.1e1 / (t32 * t28 * t20 + 0.1e1);
	t1 = [t30 * t26 * t43, 0, t22, 0, 0, 0; (-t19 * t41 - (-t24 * t28 * t30 * t42 + (t26 - 0.1e1) * t33 * t23) * t32 * t47) * t18, 0, (t35 * t19 - (-t34 * t46 + t24 * t33 + (-t24 * t41 + t46) * t22) * t47) * t36 * t18, 0, 0, 0; (-0.1e1 - t44) * t35 * t27, 0, -0.1e1 / t34 * t27 * t43, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:17:00
	% EndTime: 2019-10-10 00:17:00
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (325->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
	t50 = qJ(3) + pkin(10);
	t46 = sin(t50);
	t51 = qJ(1) + pkin(9);
	t47 = sin(t51);
	t48 = cos(t50);
	t60 = t47 * t48;
	t41 = atan2(-t60, t46);
	t39 = sin(t41);
	t40 = cos(t41);
	t33 = -t39 * t60 + t40 * t46;
	t31 = 0.1e1 / t33 ^ 2;
	t49 = cos(t51);
	t65 = t31 * t49 ^ 2;
	t52 = sin(qJ(6));
	t56 = t49 * t52;
	t53 = cos(qJ(6));
	t58 = t47 * t53;
	t38 = t46 * t56 + t58;
	t36 = 0.1e1 / t38 ^ 2;
	t55 = t49 * t53;
	t59 = t47 * t52;
	t37 = -t46 * t55 + t59;
	t64 = t36 * t37;
	t63 = t39 * t46;
	t45 = t48 ^ 2;
	t62 = 0.1e1 / t46 ^ 2 * t45;
	t42 = 0.1e1 / (t47 ^ 2 * t62 + 0.1e1);
	t61 = t47 * t42;
	t57 = t48 * t49;
	t54 = t37 ^ 2 * t36 + 0.1e1;
	t43 = 0.1e1 / t46;
	t35 = 0.1e1 / t38;
	t34 = 0.1e1 / t54;
	t32 = (0.1e1 + t62) * t61;
	t30 = 0.1e1 / t33;
	t29 = 0.1e1 / (t45 * t65 + 0.1e1);
	t1 = [-t43 * t42 * t57, 0, t32, 0, 0, 0; (-t30 * t60 - (t40 * t43 * t45 * t61 + (t42 - 0.1e1) * t48 * t39) * t48 * t65) * t29, 0, (-t46 * t30 - (t47 * t63 + t40 * t48 + (-t40 * t60 - t63) * t32) * t48 * t31) * t49 * t29, 0, 0, 0; ((t46 * t58 + t56) * t35 - (-t46 * t59 + t55) * t64) * t34, 0, (-t35 * t53 - t52 * t64) * t34 * t57, 0, 0, t54 * t34;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end