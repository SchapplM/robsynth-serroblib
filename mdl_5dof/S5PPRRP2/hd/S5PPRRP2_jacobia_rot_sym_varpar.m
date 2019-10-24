% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRRP2
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
%   Wie in S5PPRRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:19
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PPRRP2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:40
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:41
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:40
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:40
	% EndTime: 2019-10-24 10:19:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:19:40
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (137->15), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->24)
	t31 = pkin(8) + qJ(3);
	t30 = cos(t31);
	t29 = sin(t31);
	t32 = sin(pkin(7));
	t37 = t32 * t29;
	t25 = atan2(-t37, -t30);
	t23 = sin(t25);
	t40 = t23 * t30;
	t27 = t29 ^ 2;
	t39 = t27 / t30 ^ 2;
	t33 = cos(pkin(7));
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
	% StartTime: 2019-10-24 10:19:41
	% EndTime: 2019-10-24 10:19:41
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (250->22), mult. (359->55), div. (75->11), fcn. (540->9), ass. (0->37)
	t45 = pkin(8) + qJ(3);
	t42 = sin(t45);
	t62 = t42 ^ 2;
	t43 = cos(t45);
	t49 = cos(pkin(7));
	t51 = cos(qJ(4));
	t53 = t49 * t51;
	t48 = sin(pkin(7));
	t50 = sin(qJ(4));
	t56 = t48 * t50;
	t34 = t43 * t56 + t53;
	t57 = t42 * t50;
	t32 = atan2(-t34, t57);
	t28 = sin(t32);
	t29 = cos(t32);
	t27 = -t28 * t34 + t29 * t57;
	t26 = 0.1e1 / t27 ^ 2;
	t54 = t49 * t50;
	t55 = t48 * t51;
	t37 = t43 * t54 - t55;
	t61 = t26 * t37;
	t59 = t29 * t34;
	t58 = t42 * t49;
	t38 = t43 * t53 + t56;
	t33 = 0.1e1 / t38 ^ 2;
	t52 = t49 ^ 2 * t62 * t33;
	t47 = 0.1e1 / t50 ^ 2;
	t46 = 0.1e1 / t50;
	t41 = 0.1e1 / t62;
	t36 = t43 * t55 - t54;
	t31 = 0.1e1 / (t34 ^ 2 * t41 * t47 + 0.1e1);
	t30 = 0.1e1 / (0.1e1 + t52);
	t25 = 0.1e1 / t27;
	t24 = (t34 * t41 * t43 * t46 + t48) * t31;
	t23 = 0.1e1 / (t37 ^ 2 * t26 + 0.1e1);
	t22 = (t34 * t47 * t51 - t36 * t46) / t42 * t31;
	t1 = [0, 0, t24, t22, 0; 0, 0, (t24 * t59 * t61 + (-t25 * t58 - (t29 * t43 + (-t24 + t48) * t42 * t28) * t61) * t50) * t23, (t38 * t25 - (t29 * t42 * t51 - t28 * t36 + (-t28 * t57 - t59) * t22) * t61) * t23, 0; 0, 0, (-t49 * t43 / t38 - t51 * t52) * t30, -t37 * t33 * t30 * t58, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end