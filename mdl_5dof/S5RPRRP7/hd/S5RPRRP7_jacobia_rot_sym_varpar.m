% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRP7
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
%   Wie in S5RPRRP7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPRRP7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:46:17
	% EndTime: 2019-12-31 18:46:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:46:17
	% EndTime: 2019-12-31 18:46:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:46:17
	% EndTime: 2019-12-31 18:46:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:46:17
	% EndTime: 2019-12-31 18:46:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:46:17
	% EndTime: 2019-12-31 18:46:17
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (218->21), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
	t36 = qJ(1) + pkin(8);
	t35 = cos(t36);
	t54 = t35 ^ 2;
	t43 = cos(qJ(3));
	t34 = sin(t36);
	t41 = sin(qJ(3));
	t49 = t34 * t41;
	t32 = atan2(-t49, -t43);
	t30 = sin(t32);
	t31 = cos(t32);
	t23 = -t30 * t49 - t31 * t43;
	t22 = 0.1e1 / t23 ^ 2;
	t53 = t22 * t41;
	t40 = sin(qJ(4));
	t42 = cos(qJ(4));
	t45 = t42 * t43;
	t29 = t34 * t40 + t35 * t45;
	t27 = 0.1e1 / t29 ^ 2;
	t46 = t40 * t43;
	t28 = -t34 * t42 + t35 * t46;
	t52 = t27 * t28;
	t51 = t30 * t43;
	t37 = t41 ^ 2;
	t47 = t37 / t43 ^ 2;
	t33 = 0.1e1 / (t34 ^ 2 * t47 + 0.1e1);
	t50 = t34 * t33;
	t48 = t35 * t41;
	t44 = t28 ^ 2 * t27 + 0.1e1;
	t38 = 0.1e1 / t43;
	t26 = 0.1e1 / t29;
	t25 = (0.1e1 + t47) * t50;
	t24 = 0.1e1 / t44;
	t21 = 0.1e1 / t23;
	t20 = 0.1e1 / (t54 * t37 * t22 + 0.1e1);
	t1 = [t38 * t33 * t48, 0, t25, 0, 0; (-t21 * t49 - (-t31 * t37 * t38 * t50 + (t33 - 0.1e1) * t41 * t30) * t54 * t53) * t20, 0, (t43 * t21 - (-t34 * t51 + t31 * t41 + (-t31 * t49 + t51) * t25) * t53) * t35 * t20, 0, 0; ((-t34 * t46 - t35 * t42) * t26 - (-t34 * t45 + t35 * t40) * t52) * t24, 0, (-t26 * t40 + t42 * t52) * t24 * t48, t44 * t24, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:46:17
	% EndTime: 2019-12-31 18:46:17
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (378->27), mult. (509->71), div. (100->11), fcn. (770->9), ass. (0->39)
	t54 = sin(qJ(3));
	t67 = t54 ^ 2;
	t47 = qJ(1) + pkin(8);
	t45 = sin(t47);
	t46 = cos(t47);
	t55 = cos(qJ(4));
	t53 = sin(qJ(4));
	t56 = cos(qJ(3));
	t60 = t53 * t56;
	t37 = t45 * t60 + t46 * t55;
	t59 = t54 * t53;
	t34 = atan2(-t37, t59);
	t31 = sin(t34);
	t32 = cos(t34);
	t29 = -t31 * t37 + t32 * t59;
	t28 = 0.1e1 / t29 ^ 2;
	t40 = -t45 * t55 + t46 * t60;
	t66 = t28 * t40;
	t64 = t32 * t37;
	t63 = t40 ^ 2 * t28;
	t62 = t46 * t54;
	t48 = 0.1e1 / t53;
	t51 = 0.1e1 / t54;
	t61 = t48 * t51;
	t58 = t55 * t56;
	t41 = t45 * t53 + t46 * t58;
	t36 = 0.1e1 / t41 ^ 2;
	t57 = t46 ^ 2 * t67 * t36;
	t52 = 0.1e1 / t67;
	t49 = 0.1e1 / t53 ^ 2;
	t39 = t45 * t58 - t46 * t53;
	t35 = 0.1e1 / t41;
	t33 = 0.1e1 / (t37 ^ 2 * t52 * t49 + 0.1e1);
	t30 = 0.1e1 / (0.1e1 + t57);
	t27 = 0.1e1 / t29;
	t26 = (t37 * t48 * t52 * t56 + t45) * t33;
	t25 = 0.1e1 / (0.1e1 + t63);
	t24 = (t37 * t49 * t55 - t39 * t48) * t51 * t33;
	t1 = [-t40 * t33 * t61, 0, t26, t24, 0; (-t37 * t27 - (-t31 + (t61 * t64 + t31) * t33) * t63) * t25, 0, (t26 * t64 * t66 + (-t27 * t62 - (t32 * t56 + (-t26 + t45) * t54 * t31) * t66) * t53) * t25, (t41 * t27 - (t32 * t54 * t55 - t31 * t39 + (-t31 * t59 - t64) * t24) * t66) * t25, 0; (-t36 * t39 * t46 + t35 * t45) * t54 * t30, 0, (-t35 * t46 * t56 - t55 * t57) * t30, -t40 * t36 * t30 * t62, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end