% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRPRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRPRP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:39:33
	% EndTime: 2019-12-29 18:39:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:39:33
	% EndTime: 2019-12-29 18:39:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:39:38
	% EndTime: 2019-12-29 18:39:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->0), mult. (12->0), div. (10->0), fcn. (12->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:39:39
	% EndTime: 2019-12-29 18:39:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:39:33
	% EndTime: 2019-12-29 18:39:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:39:39
	% EndTime: 2019-12-29 18:39:39
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (384->18), mult. (208->39), div. (72->9), fcn. (333->7), ass. (0->31)
	t49 = qJ(1) + qJ(2);
	t47 = cos(t49);
	t50 = t47 ^ 2;
	t48 = pkin(8) + qJ(4);
	t45 = cos(t48);
	t44 = sin(t48);
	t46 = sin(t49);
	t51 = t46 * t44;
	t36 = atan2(-t51, -t45);
	t34 = sin(t36);
	t35 = cos(t36);
	t30 = -t34 * t51 - t35 * t45;
	t29 = 0.1e1 / t30 ^ 2;
	t56 = t29 * t44;
	t55 = t34 * t45;
	t39 = t44 ^ 2;
	t41 = 0.1e1 / t45 ^ 2;
	t54 = t39 * t41;
	t42 = t46 ^ 2;
	t53 = t42 / t50;
	t37 = 0.1e1 / (t42 * t54 + 0.1e1);
	t52 = t46 * t37;
	t40 = 0.1e1 / t45;
	t38 = 0.1e1 / (t41 * t53 + 0.1e1);
	t33 = t47 * t44 * t40 * t37;
	t32 = (-0.1e1 - t53) * t40 * t38;
	t31 = (0.1e1 + t54) * t52;
	t28 = 0.1e1 / t30;
	t27 = 0.1e1 / (t50 * t39 * t29 + 0.1e1);
	t26 = (-t28 * t51 - (-t35 * t39 * t40 * t52 + (t37 - 0.1e1) * t44 * t34) * t50 * t56) * t27;
	t1 = [t33, t33, 0, t31, 0; t26, t26, 0, (t45 * t28 - (-t46 * t55 + t35 * t44 + (-t35 * t51 + t55) * t31) * t56) * t47 * t27, 0; t32, t32, 0, -0.1e1 / t47 * t41 * t38 * t51, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end