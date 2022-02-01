% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPPR2
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
%   Wie in S5RPPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPPPR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
Ja_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (46->14), mult. (89->32), div. (20->9), fcn. (140->9), ass. (0->23)
	t28 = cos(pkin(7));
	t26 = sin(pkin(7));
	t29 = sin(qJ(1));
	t35 = t26 * t29;
	t21 = atan2(-t35, -t28);
	t19 = sin(t21);
	t20 = cos(t21);
	t14 = -t19 * t35 - t20 * t28;
	t30 = cos(qJ(1));
	t37 = 0.1e1 / t14 ^ 2 * t30 ^ 2;
	t23 = t26 ^ 2;
	t22 = 0.1e1 / (0.1e1 + t23 * t29 ^ 2 / t28 ^ 2);
	t36 = t22 / t28;
	t25 = sin(pkin(8));
	t34 = t29 * t25;
	t27 = cos(pkin(8));
	t33 = t29 * t27;
	t32 = t30 * t25;
	t31 = t30 * t27;
	t18 = t28 * t31 + t34;
	t17 = t28 * t32 - t33;
	t16 = 0.1e1 / t18 ^ 2;
	t1 = [t26 * t30 * t36, 0, 0, 0, 0; (-0.1e1 / t14 * t35 - (-t20 * t23 * t29 * t36 + (t22 - 0.1e1) * t26 * t19) * t26 * t37) / (t23 * t37 + 0.1e1), 0, 0, 0, 0; ((-t28 * t34 - t31) / t18 - (-t28 * t33 + t32) * t17 * t16) / (t17 ^ 2 * t16 + 0.1e1), 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (60->18), mult. (183->37), div. (25->10), fcn. (278->11), ass. (0->30)
	t41 = sin(pkin(8));
	t42 = sin(pkin(7));
	t54 = t42 * t41;
	t45 = cos(pkin(7));
	t44 = cos(pkin(8));
	t47 = cos(qJ(1));
	t48 = t47 * t44;
	t46 = sin(qJ(1));
	t51 = t46 * t41;
	t32 = t45 * t51 + t48;
	t31 = atan2(-t32, t54);
	t28 = sin(t31);
	t29 = cos(t31);
	t23 = -t28 * t32 + t29 * t54;
	t49 = t47 * t41;
	t50 = t46 * t44;
	t35 = t45 * t49 - t50;
	t56 = t35 ^ 2 / t23 ^ 2;
	t55 = 0.1e1 / t54;
	t53 = t42 * t46;
	t52 = t42 * t47;
	t43 = cos(pkin(9));
	t40 = sin(pkin(9));
	t36 = t45 * t48 + t51;
	t34 = -t45 * t50 + t49;
	t30 = 0.1e1 / (0.1e1 + t32 ^ 2 / t42 ^ 2 / t41 ^ 2);
	t27 = t36 * t43 + t40 * t52;
	t26 = t36 * t40 - t43 * t52;
	t25 = 0.1e1 / t27 ^ 2;
	t1 = [-t35 * t30 * t55, 0, 0, 0, 0; (-t32 / t23 - (-t28 + (t29 * t32 * t55 + t28) * t30) * t56) / (0.1e1 + t56), 0, 0, 0, 0; ((t34 * t40 + t43 * t53) / t27 - (t34 * t43 - t40 * t53) * t26 * t25) / (t26 ^ 2 * t25 + 0.1e1), 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (164->25), mult. (456->54), div. (26->9), fcn. (653->13), ass. (0->38)
	t68 = cos(pkin(8));
	t65 = sin(pkin(9));
	t76 = sin(pkin(7));
	t74 = t76 * t65;
	t67 = cos(pkin(9));
	t69 = cos(pkin(7));
	t78 = t69 * t67;
	t58 = t68 * t78 + t74;
	t70 = sin(qJ(5));
	t66 = sin(pkin(8));
	t72 = cos(qJ(5));
	t80 = t66 * t72;
	t51 = -t58 * t70 + t69 * t80;
	t81 = t66 * t70;
	t59 = t67 * t81 + t72 * t68;
	t71 = sin(qJ(1));
	t73 = cos(qJ(1));
	t45 = -t51 * t73 + t71 * t59;
	t50 = t58 * t72 + t69 * t81;
	t60 = t67 * t80 - t70 * t68;
	t44 = t50 * t73 + t71 * t60;
	t43 = 0.1e1 / t44 ^ 2;
	t84 = t43 * t45;
	t75 = t67 * t76;
	t79 = t68 * t69;
	t52 = (-t73 * t66 + t71 * t79) * t65 - t71 * t75;
	t57 = t68 * t74 + t78;
	t49 = atan2(-t52, t57);
	t46 = sin(t49);
	t47 = cos(t49);
	t40 = -t46 * t52 + t47 * t57;
	t54 = (t71 * t66 + t73 * t79) * t65 - t73 * t75;
	t82 = t54 ^ 2 / t40 ^ 2;
	t56 = 0.1e1 / t57;
	t48 = 0.1e1 / (0.1e1 + t52 ^ 2 / t57 ^ 2);
	t42 = 0.1e1 / t44;
	t41 = 0.1e1 / (t45 ^ 2 * t43 + 0.1e1);
	t1 = [-t54 * t56 * t48, 0, 0, 0, 0; (-t52 / t40 - (-t46 + (t47 * t52 * t56 + t46) * t48) * t82) / (0.1e1 + t82), 0, 0, 0, 0; ((t51 * t71 + t73 * t59) * t42 - (-t50 * t71 + t73 * t60) * t84) * t41, 0, 0, 0, (t44 * t42 + t45 * t84) * t41;];
	Ja_rot = t1;
end