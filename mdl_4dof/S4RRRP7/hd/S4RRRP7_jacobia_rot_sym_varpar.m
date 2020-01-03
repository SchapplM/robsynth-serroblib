% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RRRP7
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
%   Wie in S4RRRP7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S4RRRP7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRRP7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_jacobia_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:21:24
	% EndTime: 2019-12-31 17:21:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:21:24
	% EndTime: 2019-12-31 17:21:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:21:24
	% EndTime: 2019-12-31 17:21:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:21:24
	% EndTime: 2019-12-31 17:21:24
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (109->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->33)
	t40 = cos(qJ(2));
	t37 = sin(qJ(2));
	t38 = sin(qJ(1));
	t46 = t38 * t37;
	t31 = atan2(-t46, -t40);
	t29 = sin(t31);
	t30 = cos(t31);
	t22 = -t29 * t46 - t30 * t40;
	t21 = 0.1e1 / t22 ^ 2;
	t41 = cos(qJ(1));
	t51 = t21 * t41 ^ 2;
	t36 = sin(qJ(3));
	t39 = cos(qJ(3));
	t43 = t41 * t39;
	t28 = t38 * t36 + t40 * t43;
	t26 = 0.1e1 / t28 ^ 2;
	t44 = t41 * t36;
	t27 = -t38 * t39 + t40 * t44;
	t50 = t26 * t27;
	t33 = t37 ^ 2;
	t49 = t33 / t40 ^ 2;
	t48 = t37 * t41;
	t32 = 0.1e1 / (t38 ^ 2 * t49 + 0.1e1);
	t47 = t38 * t32;
	t45 = t38 * t40;
	t42 = t27 ^ 2 * t26 + 0.1e1;
	t34 = 0.1e1 / t40;
	t25 = 0.1e1 / t28;
	t24 = (0.1e1 + t49) * t47;
	t23 = 0.1e1 / t42;
	t20 = 0.1e1 / t22;
	t19 = 0.1e1 / (t33 * t51 + 0.1e1);
	t1 = [t34 * t32 * t48, t24, 0, 0; (-t20 * t46 - (-t30 * t33 * t34 * t47 + (t32 - 0.1e1) * t37 * t29) * t37 * t51) * t19, (t40 * t20 - (-t29 * t45 + t30 * t37 + (t29 * t40 - t30 * t46) * t24) * t37 * t21) * t41 * t19, 0, 0; ((-t36 * t45 - t43) * t25 - (-t39 * t45 + t44) * t50) * t23, (-t25 * t36 + t39 * t50) * t23 * t48, t42 * t23, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:21:24
	% EndTime: 2019-12-31 17:21:24
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (159->26), mult. (509->70), div. (100->11), fcn. (770->9), ass. (0->39)
	t51 = sin(qJ(2));
	t67 = t51 ^ 2;
	t50 = sin(qJ(3));
	t53 = cos(qJ(3));
	t55 = cos(qJ(1));
	t57 = t55 * t53;
	t52 = sin(qJ(1));
	t54 = cos(qJ(2));
	t59 = t52 * t54;
	t37 = t50 * t59 + t57;
	t61 = t51 * t50;
	t34 = atan2(-t37, t61);
	t30 = sin(t34);
	t31 = cos(t34);
	t29 = -t30 * t37 + t31 * t61;
	t28 = 0.1e1 / t29 ^ 2;
	t58 = t55 * t50;
	t40 = -t52 * t53 + t54 * t58;
	t66 = t28 * t40;
	t64 = t31 * t37;
	t63 = t40 ^ 2 * t28;
	t44 = 0.1e1 / t50;
	t47 = 0.1e1 / t51;
	t62 = t44 * t47;
	t60 = t51 * t55;
	t41 = t52 * t50 + t54 * t57;
	t36 = 0.1e1 / t41 ^ 2;
	t56 = t55 ^ 2 * t67 * t36;
	t48 = 0.1e1 / t67;
	t45 = 0.1e1 / t50 ^ 2;
	t39 = t53 * t59 - t58;
	t35 = 0.1e1 / t41;
	t33 = 0.1e1 / (t37 ^ 2 * t48 * t45 + 0.1e1);
	t32 = 0.1e1 / (0.1e1 + t56);
	t27 = 0.1e1 / t29;
	t26 = (t37 * t44 * t48 * t54 + t52) * t33;
	t25 = 0.1e1 / (0.1e1 + t63);
	t24 = (t37 * t45 * t53 - t39 * t44) * t47 * t33;
	t1 = [-t40 * t33 * t62, t26, t24, 0; (-t37 * t27 - (-t30 + (t62 * t64 + t30) * t33) * t63) * t25, (t26 * t64 * t66 + (-t27 * t60 - (t31 * t54 + (-t26 + t52) * t51 * t30) * t66) * t50) * t25, (t41 * t27 - (t31 * t51 * t53 - t30 * t39 + (-t30 * t61 - t64) * t24) * t66) * t25, 0; (-t36 * t39 * t55 + t35 * t52) * t51 * t32, (-t35 * t54 * t55 - t53 * t56) * t32, -t40 * t36 * t32 * t60, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,4);
end