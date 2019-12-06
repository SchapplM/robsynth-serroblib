% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPRP3
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
%   Wie in S5PRPRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRPRP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:34:00
	% EndTime: 2019-12-05 15:34:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:34:01
	% EndTime: 2019-12-05 15:34:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:34:01
	% EndTime: 2019-12-05 15:34:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:34:01
	% EndTime: 2019-12-05 15:34:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:34:01
	% EndTime: 2019-12-05 15:34:01
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (137->15), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->24)
	t32 = qJ(2) + pkin(8);
	t31 = cos(t32);
	t30 = sin(t32);
	t33 = sin(pkin(7));
	t38 = t33 * t30;
	t26 = atan2(-t38, -t31);
	t24 = sin(t26);
	t41 = t24 * t31;
	t28 = t30 ^ 2;
	t40 = t28 / t31 ^ 2;
	t34 = cos(pkin(7));
	t39 = t31 * t34;
	t35 = sin(qJ(4));
	t36 = cos(qJ(4));
	t23 = t33 * t35 + t36 * t39;
	t21 = 0.1e1 / t23 ^ 2;
	t22 = -t33 * t36 + t35 * t39;
	t37 = t22 ^ 2 * t21 + 0.1e1;
	t25 = cos(t26);
	t20 = 0.1e1 / t37;
	t19 = (0.1e1 + t40) * t33 / (t33 ^ 2 * t40 + 0.1e1);
	t18 = -t24 * t38 - t25 * t31;
	t17 = 0.1e1 / t18 ^ 2;
	t1 = [0, t19, 0, 0, 0; 0, (t31 / t18 - (-t33 * t41 + t25 * t30 + (-t25 * t38 + t41) * t19) * t30 * t17) * t34 / (t34 ^ 2 * t28 * t17 + 0.1e1), 0, 0, 0; 0, (-t35 / t23 + t36 * t22 * t21) * t34 * t30 * t20, 0, t37 * t20, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:34:01
	% EndTime: 2019-12-05 15:34:01
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (137->15), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->24)
	t37 = qJ(2) + pkin(8);
	t36 = cos(t37);
	t35 = sin(t37);
	t38 = sin(pkin(7));
	t43 = t38 * t35;
	t31 = atan2(-t43, -t36);
	t29 = sin(t31);
	t46 = t29 * t36;
	t33 = t35 ^ 2;
	t45 = t33 / t36 ^ 2;
	t39 = cos(pkin(7));
	t44 = t36 * t39;
	t40 = sin(qJ(4));
	t41 = cos(qJ(4));
	t28 = t38 * t40 + t41 * t44;
	t26 = 0.1e1 / t28 ^ 2;
	t27 = -t38 * t41 + t40 * t44;
	t42 = t27 ^ 2 * t26 + 0.1e1;
	t30 = cos(t31);
	t25 = 0.1e1 / t42;
	t24 = (0.1e1 + t45) * t38 / (t38 ^ 2 * t45 + 0.1e1);
	t23 = -t29 * t43 - t30 * t36;
	t22 = 0.1e1 / t23 ^ 2;
	t1 = [0, t24, 0, 0, 0; 0, (t36 / t23 - (-t38 * t46 + t30 * t35 + (-t30 * t43 + t46) * t24) * t35 * t22) * t39 / (t39 ^ 2 * t33 * t22 + 0.1e1), 0, 0, 0; 0, (-t40 / t28 + t41 * t27 * t26) * t39 * t35 * t25, 0, t42 * t25, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end