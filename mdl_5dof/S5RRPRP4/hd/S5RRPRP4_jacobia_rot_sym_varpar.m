% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRP4
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
%   Wie in S5RRPRP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRPRP4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_jacobia_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:53:11
	% EndTime: 2019-12-31 19:53:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:53:11
	% EndTime: 2019-12-31 19:53:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:53:11
	% EndTime: 2019-12-31 19:53:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (18->0), mult. (12->0), div. (10->0), fcn. (12->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:53:11
	% EndTime: 2019-12-31 19:53:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:53:12
	% EndTime: 2019-12-31 19:53:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:53:11
	% EndTime: 2019-12-31 19:53:11
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (188->15), mult. (208->39), div. (72->9), fcn. (333->7), ass. (0->30)
	t48 = qJ(1) + qJ(2);
	t43 = sin(t48);
	t51 = t43 ^ 2;
	t49 = sin(qJ(4));
	t44 = cos(t48);
	t50 = cos(qJ(4));
	t53 = t44 * t50;
	t39 = atan2(-t53, t49);
	t37 = sin(t39);
	t38 = cos(t39);
	t33 = -t37 * t53 + t38 * t49;
	t32 = 0.1e1 / t33 ^ 2;
	t57 = t32 * t50;
	t56 = t37 * t49;
	t42 = t44 ^ 2;
	t55 = 0.1e1 / t51 * t42;
	t46 = 0.1e1 / t49 ^ 2;
	t47 = t50 ^ 2;
	t52 = t46 * t47;
	t40 = 0.1e1 / (t42 * t52 + 0.1e1);
	t54 = t44 * t40;
	t45 = 0.1e1 / t49;
	t36 = 0.1e1 / (t46 * t55 + 0.1e1);
	t35 = t43 * t50 * t45 * t40;
	t34 = (0.1e1 + t52) * t54;
	t31 = 0.1e1 / t33;
	t30 = (0.1e1 + t55) * t45 * t36;
	t29 = 0.1e1 / (t51 * t47 * t32 + 0.1e1);
	t28 = (-t31 * t53 + (-t38 * t45 * t47 * t54 + (-t40 + 0.1e1) * t50 * t37) * t51 * t57) * t29;
	t1 = [t35, t35, 0, t34, 0; t28, t28, 0, (t49 * t31 + (t44 * t56 + t38 * t50 + (-t38 * t53 - t56) * t34) * t57) * t43 * t29, 0; t30, t30, 0, 0.1e1 / t43 * t46 * t36 * t53, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end