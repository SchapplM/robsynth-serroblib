% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S4RRRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 14:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S4RRRP6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRP6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_jacobia_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 14:02:09
	% EndTime: 2019-12-29 14:02:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 14:02:14
	% EndTime: 2019-12-29 14:02:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 14:02:14
	% EndTime: 2019-12-29 14:02:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 14:02:14
	% EndTime: 2019-12-29 14:02:14
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (109->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->33)
	t40 = cos(qJ(2));
	t37 = sin(qJ(2));
	t38 = sin(qJ(1));
	t46 = t38 * t37;
	t31 = atan2(-t46, -t40);
	t29 = sin(t31);
	t30 = cos(t31);
	t22 = -t29 * t46 - t30 * t40;
	t21 = 0.1e1 / t22 ^ 2;
	t41 = cos(qJ(1));
	t51 = t21 * t41 ^ 2;
	t36 = sin(qJ(3));
	t39 = cos(qJ(3));
	t43 = t41 * t39;
	t28 = t38 * t36 + t40 * t43;
	t26 = 0.1e1 / t28 ^ 2;
	t44 = t41 * t36;
	t27 = -t38 * t39 + t40 * t44;
	t50 = t26 * t27;
	t33 = t37 ^ 2;
	t49 = t33 / t40 ^ 2;
	t48 = t37 * t41;
	t32 = 0.1e1 / (t38 ^ 2 * t49 + 0.1e1);
	t47 = t38 * t32;
	t45 = t38 * t40;
	t42 = t27 ^ 2 * t26 + 0.1e1;
	t34 = 0.1e1 / t40;
	t25 = 0.1e1 / t28;
	t24 = (0.1e1 + t49) * t47;
	t23 = 0.1e1 / t42;
	t20 = 0.1e1 / t22;
	t19 = 0.1e1 / (t33 * t51 + 0.1e1);
	t1 = [t34 * t32 * t48, t24, 0, 0; (-t20 * t46 - (-t30 * t33 * t34 * t47 + (t32 - 0.1e1) * t37 * t29) * t37 * t51) * t19, (t40 * t20 - (-t29 * t45 + t30 * t37 + (t29 * t40 - t30 * t46) * t24) * t37 * t21) * t41 * t19, 0, 0; ((-t36 * t45 - t43) * t25 - (-t39 * t45 + t44) * t50) * t23, (-t25 * t36 + t39 * t50) * t23 * t48, t42 * t23, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 14:02:14
	% EndTime: 2019-12-29 14:02:14
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (109->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->33)
	t45 = cos(qJ(2));
	t42 = sin(qJ(2));
	t43 = sin(qJ(1));
	t51 = t43 * t42;
	t36 = atan2(-t51, -t45);
	t34 = sin(t36);
	t35 = cos(t36);
	t27 = -t34 * t51 - t35 * t45;
	t26 = 0.1e1 / t27 ^ 2;
	t46 = cos(qJ(1));
	t56 = t26 * t46 ^ 2;
	t41 = sin(qJ(3));
	t44 = cos(qJ(3));
	t48 = t46 * t44;
	t33 = t43 * t41 + t45 * t48;
	t31 = 0.1e1 / t33 ^ 2;
	t49 = t46 * t41;
	t32 = -t43 * t44 + t45 * t49;
	t55 = t31 * t32;
	t38 = t42 ^ 2;
	t54 = t38 / t45 ^ 2;
	t53 = t42 * t46;
	t37 = 0.1e1 / (t43 ^ 2 * t54 + 0.1e1);
	t52 = t43 * t37;
	t50 = t43 * t45;
	t47 = t32 ^ 2 * t31 + 0.1e1;
	t39 = 0.1e1 / t45;
	t30 = 0.1e1 / t33;
	t29 = (0.1e1 + t54) * t52;
	t28 = 0.1e1 / t47;
	t25 = 0.1e1 / t27;
	t24 = 0.1e1 / (t38 * t56 + 0.1e1);
	t1 = [t39 * t37 * t53, t29, 0, 0; (-t25 * t51 - (-t35 * t38 * t39 * t52 + (t37 - 0.1e1) * t42 * t34) * t42 * t56) * t24, (t45 * t25 - (-t34 * t50 + t35 * t42 + (t34 * t45 - t35 * t51) * t29) * t42 * t26) * t46 * t24, 0, 0; ((-t41 * t50 - t48) * t30 - (-t44 * t50 + t49) * t55) * t28, (-t30 * t41 + t44 * t55) * t28 * t53, t47 * t28, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,4);
end