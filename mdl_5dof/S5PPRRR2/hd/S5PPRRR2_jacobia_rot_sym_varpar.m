% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRRR2
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
%   Wie in S5PPRRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:20
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PPRRR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:55
	% EndTime: 2019-10-24 10:20:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:55
	% EndTime: 2019-10-24 10:20:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:55
	% EndTime: 2019-10-24 10:20:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:55
	% EndTime: 2019-10-24 10:20:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:55
	% EndTime: 2019-10-24 10:20:55
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (137->15), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->24)
	t31 = pkin(9) + qJ(3);
	t30 = cos(t31);
	t29 = sin(t31);
	t32 = sin(pkin(8));
	t37 = t32 * t29;
	t25 = atan2(-t37, -t30);
	t23 = sin(t25);
	t40 = t23 * t30;
	t27 = t29 ^ 2;
	t39 = t27 / t30 ^ 2;
	t33 = cos(pkin(8));
	t38 = t30 * t33;
	t34 = sin(qJ(4));
	t35 = cos(qJ(4));
	t22 = t32 * t34 + t35 * t38;
	t20 = 0.1e1 / t22 ^ 2;
	t21 = -t32 * t35 + t34 * t38;
	t36 = t21 ^ 2 * t20 + 0.1e1;
	t24 = cos(t25);
	t19 = 0.1e1 / t36;
	t18 = (0.1e1 + t39) * t32 / (t32 ^ 2 * t39 + 0.1e1);
	t17 = -t23 * t37 - t24 * t30;
	t16 = 0.1e1 / t17 ^ 2;
	t1 = [0, 0, t18, 0, 0; 0, 0, (t30 / t17 - (-t32 * t40 + t24 * t29 + (-t24 * t37 + t40) * t18) * t29 * t16) * t33 / (t33 ^ 2 * t27 * t16 + 0.1e1), 0, 0; 0, 0, (-t34 / t22 + t35 * t21 * t20) * t33 * t29 * t19, t36 * t19, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:55
	% EndTime: 2019-10-24 10:20:55
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (198->16), mult. (162->34), div. (37->8), fcn. (227->9), ass. (0->26)
	t47 = pkin(9) + qJ(3);
	t44 = cos(t47);
	t43 = sin(t47);
	t49 = sin(pkin(8));
	t52 = t49 * t43;
	t39 = atan2(-t52, -t44);
	t37 = sin(t39);
	t55 = t37 * t44;
	t41 = t43 ^ 2;
	t54 = t41 / t44 ^ 2;
	t50 = cos(pkin(8));
	t53 = t44 * t50;
	t48 = qJ(4) + qJ(5);
	t45 = sin(t48);
	t46 = cos(t48);
	t36 = t49 * t45 + t46 * t53;
	t34 = 0.1e1 / t36 ^ 2;
	t35 = t45 * t53 - t49 * t46;
	t51 = t35 ^ 2 * t34 + 0.1e1;
	t38 = cos(t39);
	t33 = 0.1e1 / t51;
	t32 = (0.1e1 + t54) * t49 / (t49 ^ 2 * t54 + 0.1e1);
	t31 = -t37 * t52 - t38 * t44;
	t30 = 0.1e1 / t31 ^ 2;
	t28 = t51 * t33;
	t1 = [0, 0, t32, 0, 0; 0, 0, (t44 / t31 - (-t49 * t55 + t38 * t43 + (-t38 * t52 + t55) * t32) * t43 * t30) * t50 / (t50 ^ 2 * t41 * t30 + 0.1e1), 0, 0; 0, 0, (-t45 / t36 + t46 * t35 * t34) * t50 * t43 * t33, t28, t28;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end