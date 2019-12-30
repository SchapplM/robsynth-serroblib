% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRPP1
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
%   Wie in S5RPRPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPRPP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:32:29
	% EndTime: 2019-12-29 16:32:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:32:29
	% EndTime: 2019-12-29 16:32:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:32:29
	% EndTime: 2019-12-29 16:32:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:32:29
	% EndTime: 2019-12-29 16:32:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:32:29
	% EndTime: 2019-12-29 16:32:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:32:29
	% EndTime: 2019-12-29 16:32:29
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (264->18), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->28)
	t38 = qJ(1) + pkin(7);
	t36 = cos(t38);
	t39 = t36 ^ 2;
	t37 = qJ(3) + pkin(8);
	t35 = cos(t37);
	t33 = sin(t37);
	t34 = sin(t38);
	t40 = t34 * t33;
	t25 = atan2(-t40, -t35);
	t23 = sin(t25);
	t24 = cos(t25);
	t21 = -t23 * t40 - t24 * t35;
	t20 = 0.1e1 / t21 ^ 2;
	t45 = t20 * t33;
	t44 = t23 * t35;
	t28 = t33 ^ 2;
	t31 = 0.1e1 / t35 ^ 2;
	t43 = t28 * t31;
	t29 = t34 ^ 2;
	t42 = t29 / t39;
	t26 = 0.1e1 / (t29 * t43 + 0.1e1);
	t41 = t34 * t26;
	t30 = 0.1e1 / t35;
	t27 = 0.1e1 / (t31 * t42 + 0.1e1);
	t22 = (0.1e1 + t43) * t41;
	t19 = 0.1e1 / t21;
	t18 = 0.1e1 / (t39 * t28 * t20 + 0.1e1);
	t1 = [t36 * t33 * t30 * t26, 0, t22, 0, 0; (-t19 * t40 - (-t24 * t28 * t30 * t41 + (t26 - 0.1e1) * t33 * t23) * t39 * t45) * t18, 0, (t35 * t19 - (-t34 * t44 + t24 * t33 + (-t24 * t40 + t44) * t22) * t45) * t36 * t18, 0, 0; (-0.1e1 - t42) * t30 * t27, 0, -0.1e1 / t36 * t31 * t27 * t40, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end