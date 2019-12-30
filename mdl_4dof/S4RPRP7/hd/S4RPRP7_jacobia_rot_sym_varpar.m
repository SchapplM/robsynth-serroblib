% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RPRP7
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
%   Wie in S4RPRP7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S4RPRP7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPRP7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_jacobia_rot_sym_varpar: pkin has to be [5x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:00:30
	% EndTime: 2019-12-29 13:00:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:00:30
	% EndTime: 2019-12-29 13:00:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:00:30
	% EndTime: 2019-12-29 13:00:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:00:30
	% EndTime: 2019-12-29 13:00:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:00:30
	% EndTime: 2019-12-29 13:00:30
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (58->14), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->26)
	t27 = sin(qJ(1));
	t30 = t27 ^ 2;
	t26 = sin(qJ(3));
	t28 = cos(qJ(3));
	t29 = cos(qJ(1));
	t31 = t29 * t28;
	t18 = atan2(-t31, t26);
	t16 = sin(t18);
	t17 = cos(t18);
	t14 = -t16 * t31 + t17 * t26;
	t13 = 0.1e1 / t14 ^ 2;
	t36 = t13 * t28;
	t35 = t16 * t26;
	t22 = 0.1e1 / t26 ^ 2;
	t24 = t28 ^ 2;
	t34 = t22 * t24;
	t25 = t29 ^ 2;
	t33 = 0.1e1 / t30 * t25;
	t20 = 0.1e1 / (t25 * t34 + 0.1e1);
	t32 = t29 * t20;
	t21 = 0.1e1 / t26;
	t19 = 0.1e1 / (t22 * t33 + 0.1e1);
	t15 = (0.1e1 + t34) * t32;
	t12 = 0.1e1 / t14;
	t11 = 0.1e1 / (t30 * t24 * t13 + 0.1e1);
	t1 = [t27 * t28 * t21 * t20, 0, t15, 0; (-t12 * t31 + (-t17 * t21 * t24 * t32 + (-t20 + 0.1e1) * t28 * t16) * t30 * t36) * t11, 0, (t26 * t12 + (t29 * t35 + t17 * t28 + (-t17 * t31 - t35) * t15) * t36) * t27 * t11, 0; (0.1e1 + t33) * t21 * t19, 0, 0.1e1 / t27 * t22 * t19 * t31, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,4);
end