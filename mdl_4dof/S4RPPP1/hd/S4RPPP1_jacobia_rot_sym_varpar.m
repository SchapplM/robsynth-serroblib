% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RPPP1
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
%   Wie in S4RPPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S4RPPP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobia_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (25->13), mult. (89->32), div. (20->9), fcn. (140->9), ass. (0->23)
	t28 = cos(pkin(4));
	t26 = sin(pkin(4));
	t30 = cos(qJ(1));
	t32 = t30 * t26;
	t22 = atan2(t32, t28);
	t19 = sin(t22);
	t20 = cos(t22);
	t14 = t19 * t32 + t20 * t28;
	t29 = sin(qJ(1));
	t37 = 0.1e1 / t14 ^ 2 * t29 ^ 2;
	t23 = t26 ^ 2;
	t21 = 0.1e1 / (0.1e1 + t30 ^ 2 * t23 / t28 ^ 2);
	t36 = t21 / t28;
	t25 = sin(pkin(6));
	t35 = t29 * t25;
	t27 = cos(pkin(6));
	t34 = t29 * t27;
	t33 = t30 * t25;
	t31 = t30 * t27;
	t18 = -t28 * t35 + t31;
	t17 = t28 * t34 + t33;
	t16 = 0.1e1 / t18 ^ 2;
	t1 = [-t29 * t26 * t36, 0, 0, 0; (0.1e1 / t14 * t32 - (-t20 * t23 * t30 * t36 + (t21 - 0.1e1) * t26 * t19) * t26 * t37) / (t23 * t37 + 0.1e1), 0, 0, 0; ((t28 * t31 - t35) / t18 - (-t28 * t33 - t34) * t17 * t16) / (t17 ^ 2 * t16 + 0.1e1), 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (58->14), mult. (143->30), div. (29->11), fcn. (226->9), ass. (0->23)
	t37 = cos(pkin(4));
	t36 = cos(pkin(6));
	t39 = cos(qJ(1));
	t40 = t39 * t36;
	t34 = sin(pkin(6));
	t38 = sin(qJ(1));
	t41 = t38 * t34;
	t25 = -t37 * t40 + t41;
	t35 = sin(pkin(4));
	t42 = t35 * t36;
	t22 = atan2(-t25, -t42);
	t20 = sin(t22);
	t21 = cos(t22);
	t19 = -t20 * t25 - t21 * t42;
	t27 = t38 * t37 * t36 + t39 * t34;
	t44 = t27 ^ 2 / t19 ^ 2;
	t30 = 0.1e1 / t35;
	t43 = t30 / t36;
	t33 = 0.1e1 / t38 ^ 2;
	t31 = 0.1e1 / t35 ^ 2;
	t28 = -t37 * t41 + t40;
	t23 = 0.1e1 / (0.1e1 + t25 ^ 2 * t31 / t36 ^ 2);
	t1 = [t27 * t23 * t43, 0, 0, 0; (-t25 / t19 - (-t20 + (-t21 * t25 * t43 + t20) * t23) * t44) / (0.1e1 + t44), 0, 0, 0; (-t36 + (-t34 * t37 / t38 - t28 * t33) * t39) * t30 / (t28 ^ 2 * t33 * t31 + 0.1e1), 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (46->14), mult. (143->30), div. (29->11), fcn. (226->9), ass. (0->23)
	t37 = cos(pkin(4));
	t34 = sin(pkin(6));
	t39 = cos(qJ(1));
	t40 = t39 * t34;
	t36 = cos(pkin(6));
	t38 = sin(qJ(1));
	t41 = t38 * t36;
	t25 = t37 * t40 + t41;
	t35 = sin(pkin(4));
	t42 = t35 * t34;
	t23 = atan2(-t25, t42);
	t20 = sin(t23);
	t21 = cos(t23);
	t19 = -t20 * t25 + t21 * t42;
	t28 = -t38 * t37 * t34 + t39 * t36;
	t44 = t28 ^ 2 / t19 ^ 2;
	t31 = 0.1e1 / t35;
	t43 = 0.1e1 / t34 * t31;
	t33 = 0.1e1 / t38 ^ 2;
	t32 = 0.1e1 / t35 ^ 2;
	t27 = -t37 * t41 - t40;
	t22 = 0.1e1 / (0.1e1 + t25 ^ 2 * t32 / t34 ^ 2);
	t1 = [-t28 * t22 * t43, 0, 0, 0; (-t25 / t19 - (-t20 + (t21 * t25 * t43 + t20) * t22) * t44) / (0.1e1 + t44), 0, 0, 0; (t34 + (-t36 * t37 / t38 - t27 * t33) * t39) * t31 / (t27 ^ 2 * t33 * t32 + 0.1e1), 0, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,4);
end