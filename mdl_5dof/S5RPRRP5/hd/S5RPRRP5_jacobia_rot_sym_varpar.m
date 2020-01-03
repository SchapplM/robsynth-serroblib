% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRP5
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
%   Wie in S5RPRRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPRRP5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:41:20
	% EndTime: 2019-12-31 18:41:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:41:20
	% EndTime: 2019-12-31 18:41:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:41:20
	% EndTime: 2019-12-31 18:41:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:41:20
	% EndTime: 2019-12-31 18:41:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (30->0), mult. (12->0), div. (10->0), fcn. (12->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 1, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:41:20
	% EndTime: 2019-12-31 18:41:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:41:20
	% EndTime: 2019-12-31 18:41:20
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (327->17), mult. (208->39), div. (72->9), fcn. (333->7), ass. (0->30)
	t43 = qJ(1) + pkin(8) + qJ(3);
	t42 = cos(t43);
	t49 = t42 ^ 2;
	t48 = cos(qJ(4));
	t41 = sin(t43);
	t47 = sin(qJ(4));
	t51 = t41 * t47;
	t37 = atan2(-t51, -t48);
	t35 = sin(t37);
	t36 = cos(t37);
	t31 = -t35 * t51 - t36 * t48;
	t30 = 0.1e1 / t31 ^ 2;
	t55 = t30 * t47;
	t54 = t35 * t48;
	t39 = t41 ^ 2;
	t53 = t39 / t49;
	t44 = t47 ^ 2;
	t46 = 0.1e1 / t48 ^ 2;
	t50 = t44 * t46;
	t38 = 0.1e1 / (t39 * t50 + 0.1e1);
	t52 = t41 * t38;
	t45 = 0.1e1 / t48;
	t34 = 0.1e1 / (t46 * t53 + 0.1e1);
	t33 = t42 * t47 * t45 * t38;
	t32 = (0.1e1 + t50) * t52;
	t29 = 0.1e1 / t31;
	t28 = (-0.1e1 - t53) * t45 * t34;
	t27 = 0.1e1 / (t49 * t44 * t30 + 0.1e1);
	t26 = (-t29 * t51 - (-t36 * t44 * t45 * t52 + (t38 - 0.1e1) * t47 * t35) * t49 * t55) * t27;
	t1 = [t33, 0, t33, t32, 0; t26, 0, t26, (t48 * t29 - (-t41 * t54 + t36 * t47 + (-t36 * t51 + t54) * t32) * t55) * t42 * t27, 0; t28, 0, t28, -0.1e1 / t42 * t46 * t34 * t51, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end