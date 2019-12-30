% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RRPP3
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
%   Wie in S4RRPP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:23
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S4RRPP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_jacobia_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:23:30
	% EndTime: 2019-12-29 13:23:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:23:37
	% EndTime: 2019-12-29 13:23:37
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:23:37
	% EndTime: 2019-12-29 13:23:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:23:37
	% EndTime: 2019-12-29 13:23:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:23:37
	% EndTime: 2019-12-29 13:23:37
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (193->17), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->27)
	t36 = cos(qJ(1));
	t37 = t36 ^ 2;
	t32 = qJ(2) + pkin(6);
	t31 = cos(t32);
	t30 = sin(t32);
	t35 = sin(qJ(1));
	t38 = t35 * t30;
	t24 = atan2(-t38, -t31);
	t22 = sin(t24);
	t23 = cos(t24);
	t20 = -t22 * t38 - t23 * t31;
	t19 = 0.1e1 / t20 ^ 2;
	t43 = t19 * t30;
	t42 = t22 * t31;
	t27 = t30 ^ 2;
	t29 = 0.1e1 / t31 ^ 2;
	t41 = t27 * t29;
	t33 = t35 ^ 2;
	t40 = t33 / t37;
	t25 = 0.1e1 / (t33 * t41 + 0.1e1);
	t39 = t35 * t25;
	t28 = 0.1e1 / t31;
	t26 = 0.1e1 / (t29 * t40 + 0.1e1);
	t21 = (0.1e1 + t41) * t39;
	t18 = 0.1e1 / t20;
	t17 = 0.1e1 / (t37 * t27 * t19 + 0.1e1);
	t1 = [t36 * t30 * t28 * t25, t21, 0, 0; (-t18 * t38 - (-t23 * t27 * t28 * t39 + (t25 - 0.1e1) * t30 * t22) * t37 * t43) * t17, (t31 * t18 - (-t35 * t42 + t23 * t30 + (-t23 * t38 + t42) * t21) * t43) * t36 * t17, 0, 0; (-0.1e1 - t40) * t28 * t26, -0.1e1 / t36 * t29 * t26 * t38, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,4);
end