% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRPR3
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
%   Wie in S5PPRPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PPRPR3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 1, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (24->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 1, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (222->18), mult. (268->42), div. (47->11), fcn. (395->11), ass. (0->31)
	t45 = qJ(3) + pkin(9);
	t42 = sin(t45);
	t46 = sin(pkin(8));
	t57 = t46 * t42;
	t49 = cos(pkin(7));
	t56 = t46 * t49;
	t47 = sin(pkin(7));
	t48 = cos(pkin(8));
	t55 = t47 * t48;
	t54 = t49 * t42;
	t43 = cos(t45);
	t53 = t49 * t43;
	t39 = t47 * t42 + t48 * t53;
	t50 = sin(qJ(5));
	t51 = cos(qJ(5));
	t30 = t39 * t51 + t50 * t56;
	t28 = 0.1e1 / t30 ^ 2;
	t29 = t39 * t50 - t51 * t56;
	t52 = t29 ^ 2 * t28 + 0.1e1;
	t41 = 0.1e1 / t42 ^ 2;
	t38 = -t47 * t43 + t48 * t54;
	t37 = t43 * t55 - t54;
	t35 = t42 * t55 + t53;
	t34 = atan2(-t35, t57);
	t32 = cos(t34);
	t31 = sin(t34);
	t27 = 0.1e1 / t52;
	t26 = -t31 * t35 + t32 * t57;
	t25 = 0.1e1 / t26 ^ 2;
	t23 = (-t37 / t42 + t43 * t35 * t41) / t46 / (0.1e1 + t35 ^ 2 / t46 ^ 2 * t41);
	t1 = [0, 0, t23, 0, 0; 0, 0, (t39 / t26 - (t32 * t46 * t43 - t31 * t37 + (-t31 * t57 - t32 * t35) * t23) * t38 * t25) / (t38 ^ 2 * t25 + 0.1e1), 0, 0; 0, 0, (-t50 / t30 + t51 * t29 * t28) * t38 * t27, 0, t52 * t27;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end