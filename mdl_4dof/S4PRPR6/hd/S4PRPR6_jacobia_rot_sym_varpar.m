% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4PRPR6
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
%   Wie in S4PRPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S4PRPR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRPR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_jacobia_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:09:29
	% EndTime: 2019-12-29 12:09:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:09:29
	% EndTime: 2019-12-29 12:09:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:09:29
	% EndTime: 2019-12-29 12:09:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:09:24
	% EndTime: 2019-12-29 12:09:24
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (53->14), mult. (108->32), div. (27->8), fcn. (157->9), ass. (0->21)
	t29 = cos(qJ(2));
	t25 = sin(pkin(6));
	t28 = sin(qJ(2));
	t31 = t25 * t28;
	t20 = atan2(-t31, -t29);
	t18 = sin(t20);
	t33 = t18 * t29;
	t22 = t28 ^ 2;
	t32 = t22 / t29 ^ 2;
	t27 = cos(pkin(6));
	t30 = t27 * t29;
	t26 = cos(pkin(7));
	t24 = sin(pkin(7));
	t19 = cos(t20);
	t17 = t25 * t24 + t26 * t30;
	t16 = t24 * t30 - t25 * t26;
	t15 = 0.1e1 / t17 ^ 2;
	t14 = (0.1e1 + t32) * t25 / (t25 ^ 2 * t32 + 0.1e1);
	t12 = -t18 * t31 - t19 * t29;
	t11 = 0.1e1 / t12 ^ 2;
	t1 = [0, t14, 0, 0; 0, (t29 / t12 - (-t25 * t33 + t19 * t28 + (-t19 * t31 + t33) * t14) * t28 * t11) * t27 / (t27 ^ 2 * t22 * t11 + 0.1e1), 0, 0; 0, (-t24 / t17 + t26 * t16 * t15) * t28 * t27 / (t16 ^ 2 * t15 + 0.1e1), 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:09:29
	% EndTime: 2019-12-29 12:09:29
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (93->15), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->24)
	t35 = cos(qJ(2));
	t32 = sin(pkin(6));
	t34 = sin(qJ(2));
	t38 = t32 * t34;
	t25 = atan2(-t38, -t35);
	t23 = sin(t25);
	t40 = t23 * t35;
	t30 = t34 ^ 2;
	t39 = t30 / t35 ^ 2;
	t33 = cos(pkin(6));
	t37 = t33 * t35;
	t29 = pkin(7) + qJ(4);
	t27 = sin(t29);
	t28 = cos(t29);
	t22 = t32 * t27 + t28 * t37;
	t20 = 0.1e1 / t22 ^ 2;
	t21 = t27 * t37 - t32 * t28;
	t36 = t21 ^ 2 * t20 + 0.1e1;
	t24 = cos(t25);
	t19 = (0.1e1 + t39) * t32 / (t32 ^ 2 * t39 + 0.1e1);
	t18 = -t23 * t38 - t24 * t35;
	t17 = 0.1e1 / t18 ^ 2;
	t16 = 0.1e1 / t36;
	t1 = [0, t19, 0, 0; 0, (t35 / t18 - (-t32 * t40 + t24 * t34 + (-t24 * t38 + t40) * t19) * t34 * t17) * t33 / (t33 ^ 2 * t30 * t17 + 0.1e1), 0, 0; 0, (-t27 / t22 + t28 * t21 * t20) * t34 * t33 * t16, 0, t36 * t16;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,4);
end