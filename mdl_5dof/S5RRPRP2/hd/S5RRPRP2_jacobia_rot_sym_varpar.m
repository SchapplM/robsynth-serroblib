% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRP2
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
%   Wie in S5RRPRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:37
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRPRP2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:36:54
	% EndTime: 2019-12-29 18:36:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:36:54
	% EndTime: 2019-12-29 18:36:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:36:55
	% EndTime: 2019-12-29 18:36:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->0), mult. (12->0), div. (10->0), fcn. (12->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:36:55
	% EndTime: 2019-12-29 18:36:55
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->0), mult. (12->0), div. (10->0), fcn. (12->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:36:55
	% EndTime: 2019-12-29 18:36:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:36:55
	% EndTime: 2019-12-29 18:36:55
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (327->17), mult. (208->39), div. (72->9), fcn. (333->7), ass. (0->30)
	t45 = qJ(1) + qJ(2) + pkin(8);
	t44 = cos(t45);
	t51 = t44 ^ 2;
	t50 = cos(qJ(4));
	t43 = sin(t45);
	t49 = sin(qJ(4));
	t53 = t43 * t49;
	t39 = atan2(-t53, -t50);
	t37 = sin(t39);
	t38 = cos(t39);
	t33 = -t37 * t53 - t38 * t50;
	t32 = 0.1e1 / t33 ^ 2;
	t57 = t32 * t49;
	t56 = t37 * t50;
	t41 = t43 ^ 2;
	t55 = t41 / t51;
	t46 = t49 ^ 2;
	t48 = 0.1e1 / t50 ^ 2;
	t52 = t46 * t48;
	t40 = 0.1e1 / (t41 * t52 + 0.1e1);
	t54 = t43 * t40;
	t47 = 0.1e1 / t50;
	t36 = 0.1e1 / (t48 * t55 + 0.1e1);
	t35 = t44 * t49 * t47 * t40;
	t34 = (0.1e1 + t52) * t54;
	t31 = 0.1e1 / t33;
	t30 = (-0.1e1 - t55) * t47 * t36;
	t29 = 0.1e1 / (t51 * t46 * t32 + 0.1e1);
	t28 = (-t31 * t53 - (-t38 * t46 * t47 * t54 + (t40 - 0.1e1) * t49 * t37) * t51 * t57) * t29;
	t1 = [t35, t35, 0, t34, 0; t28, t28, 0, (t50 * t31 - (-t43 * t56 + t38 * t49 + (-t38 * t53 + t56) * t34) * t57) * t44 * t29, 0; t30, t30, 0, -0.1e1 / t44 * t48 * t36 * t53, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end