% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPRP6
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
%   Wie in S5RPPRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPPRP6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_jacobia_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:55:42
	% EndTime: 2019-12-31 17:55:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:55:42
	% EndTime: 2019-12-31 17:55:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:55:42
	% EndTime: 2019-12-31 17:55:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:55:42
	% EndTime: 2019-12-31 17:55:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:55:42
	% EndTime: 2019-12-31 17:55:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:55:42
	% EndTime: 2019-12-31 17:55:42
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (169->15), mult. (145->41), div. (49->9), fcn. (228->7), ass. (0->27)
	t32 = sin(qJ(1));
	t34 = t32 ^ 2;
	t29 = pkin(7) + qJ(4);
	t27 = sin(t29);
	t28 = cos(t29);
	t33 = cos(qJ(1));
	t35 = t33 * t28;
	t21 = atan2(-t35, t27);
	t19 = sin(t21);
	t20 = cos(t21);
	t17 = -t19 * t35 + t20 * t27;
	t16 = 0.1e1 / t17 ^ 2;
	t40 = t16 * t28;
	t39 = t19 * t27;
	t25 = 0.1e1 / t27 ^ 2;
	t26 = t28 ^ 2;
	t38 = t25 * t26;
	t31 = t33 ^ 2;
	t37 = 0.1e1 / t34 * t31;
	t22 = 0.1e1 / (t31 * t38 + 0.1e1);
	t36 = t33 * t22;
	t24 = 0.1e1 / t27;
	t23 = 0.1e1 / (t25 * t37 + 0.1e1);
	t18 = (0.1e1 + t38) * t36;
	t15 = 0.1e1 / t17;
	t14 = 0.1e1 / (t34 * t26 * t16 + 0.1e1);
	t1 = [t32 * t28 * t24 * t22, 0, 0, t18, 0; (-t15 * t35 + (-t20 * t24 * t26 * t36 + (-t22 + 0.1e1) * t28 * t19) * t34 * t40) * t14, 0, 0, (t27 * t15 + (t33 * t39 + t20 * t28 + (-t20 * t35 - t39) * t18) * t40) * t32 * t14, 0; (0.1e1 + t37) * t24 * t23, 0, 0, 0.1e1 / t32 * t25 * t23 * t35, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end