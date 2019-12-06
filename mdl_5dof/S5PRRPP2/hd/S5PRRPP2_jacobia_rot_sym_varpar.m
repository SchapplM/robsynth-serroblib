% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRPP2
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
%   Wie in S5PRRPP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRRPP2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:55
	% EndTime: 2019-12-05 16:10:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:55
	% EndTime: 2019-12-05 16:10:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:55
	% EndTime: 2019-12-05 16:10:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:55
	% EndTime: 2019-12-05 16:10:56
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (63->14), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->23)
	t31 = cos(qJ(2));
	t26 = sin(pkin(7));
	t29 = sin(qJ(2));
	t34 = t26 * t29;
	t22 = atan2(-t34, -t31);
	t20 = sin(t22);
	t36 = t20 * t31;
	t24 = t29 ^ 2;
	t35 = t24 / t31 ^ 2;
	t27 = cos(pkin(7));
	t33 = t27 * t31;
	t28 = sin(qJ(3));
	t30 = cos(qJ(3));
	t19 = t26 * t28 + t30 * t33;
	t17 = 0.1e1 / t19 ^ 2;
	t18 = -t26 * t30 + t28 * t33;
	t32 = t18 ^ 2 * t17 + 0.1e1;
	t21 = cos(t22);
	t16 = (0.1e1 + t35) * t26 / (t26 ^ 2 * t35 + 0.1e1);
	t15 = 0.1e1 / t32;
	t14 = -t20 * t34 - t21 * t31;
	t13 = 0.1e1 / t14 ^ 2;
	t1 = [0, t16, 0, 0, 0; 0, (t31 / t14 - (-t26 * t36 + t21 * t29 + (-t21 * t34 + t36) * t16) * t29 * t13) * t27 / (t27 ^ 2 * t24 * t13 + 0.1e1), 0, 0, 0; 0, (-t28 / t19 + t30 * t18 * t17) * t29 * t27 * t15, t32 * t15, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:55
	% EndTime: 2019-12-05 16:10:55
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (93->15), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->24)
	t38 = cos(qJ(2));
	t35 = sin(pkin(7));
	t37 = sin(qJ(2));
	t41 = t35 * t37;
	t28 = atan2(-t41, -t38);
	t26 = sin(t28);
	t43 = t26 * t38;
	t33 = t37 ^ 2;
	t42 = t33 / t38 ^ 2;
	t36 = cos(pkin(7));
	t40 = t36 * t38;
	t32 = qJ(3) + pkin(8);
	t30 = sin(t32);
	t31 = cos(t32);
	t25 = t35 * t30 + t31 * t40;
	t23 = 0.1e1 / t25 ^ 2;
	t24 = t30 * t40 - t35 * t31;
	t39 = t24 ^ 2 * t23 + 0.1e1;
	t27 = cos(t28);
	t22 = (0.1e1 + t42) * t35 / (t35 ^ 2 * t42 + 0.1e1);
	t21 = -t26 * t41 - t27 * t38;
	t20 = 0.1e1 / t21 ^ 2;
	t19 = 0.1e1 / t39;
	t1 = [0, t22, 0, 0, 0; 0, (t38 / t21 - (-t35 * t43 + t27 * t37 + (-t27 * t41 + t43) * t22) * t37 * t20) * t36 / (t36 ^ 2 * t33 * t20 + 0.1e1), 0, 0, 0; 0, (-t30 / t25 + t31 * t24 * t23) * t37 * t36 * t19, t39 * t19, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:55
	% EndTime: 2019-12-05 16:10:56
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (315->22), mult. (359->56), div. (75->11), fcn. (540->9), ass. (0->35)
	t53 = sin(qJ(2));
	t63 = t53 ^ 2;
	t47 = qJ(3) + pkin(8);
	t44 = sin(t47);
	t45 = cos(t47);
	t52 = cos(pkin(7));
	t51 = sin(pkin(7));
	t54 = cos(qJ(2));
	t59 = t51 * t54;
	t37 = t44 * t59 + t52 * t45;
	t56 = t53 * t44;
	t34 = atan2(-t37, t56);
	t31 = sin(t34);
	t32 = cos(t34);
	t30 = -t31 * t37 + t32 * t56;
	t29 = 0.1e1 / t30 ^ 2;
	t57 = t52 * t54;
	t40 = t44 * t57 - t51 * t45;
	t62 = t29 * t40;
	t60 = t32 * t37;
	t58 = t52 * t53;
	t41 = t51 * t44 + t45 * t57;
	t36 = 0.1e1 / t41 ^ 2;
	t55 = t52 ^ 2 * t63 * t36;
	t50 = 0.1e1 / t63;
	t43 = 0.1e1 / t44 ^ 2;
	t42 = 0.1e1 / t44;
	t39 = -t52 * t44 + t45 * t59;
	t35 = 0.1e1 / (0.1e1 + t55);
	t33 = 0.1e1 / (t37 ^ 2 * t50 * t43 + 0.1e1);
	t28 = 0.1e1 / t30;
	t27 = (t37 * t42 * t50 * t54 + t51) * t33;
	t26 = 0.1e1 / (t40 ^ 2 * t29 + 0.1e1);
	t25 = (t37 * t43 * t45 - t39 * t42) / t53 * t33;
	t1 = [0, t27, t25, 0, 0; 0, (t27 * t60 * t62 + (-t28 * t58 - (t32 * t54 + (-t27 + t51) * t53 * t31) * t62) * t44) * t26, (t41 * t28 - (t32 * t53 * t45 - t31 * t39 + (-t31 * t56 - t60) * t25) * t62) * t26, 0, 0; 0, (-0.1e1 / t41 * t57 - t45 * t55) * t35, -t40 * t36 * t35 * t58, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end