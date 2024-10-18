% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRRR11
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
%   Wie in S5PRRRR11_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRRRR11_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR11_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
Ja_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	% Symbolic code from jacobia_rot_1_floatb_twist_matlab.m not found
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (94->14), mult. (116->36), div. (25->9), fcn. (175->9), ass. (0->24)
	t34 = cos(pkin(5));
	t32 = pkin(10) + qJ(2);
	t29 = cos(t32);
	t33 = sin(pkin(5));
	t40 = t29 * t33;
	t27 = atan2(t40, t34);
	t24 = sin(t27);
	t25 = cos(t27);
	t20 = t24 * t40 + t25 * t34;
	t28 = sin(t32);
	t42 = 0.1e1 / t20 ^ 2 * t28 ^ 2;
	t30 = t33 ^ 2;
	t26 = 0.1e1 / (0.1e1 + t29 ^ 2 * t30 / t34 ^ 2);
	t41 = t26 / t34;
	t35 = sin(qJ(3));
	t39 = t34 * t35;
	t36 = cos(qJ(3));
	t38 = t34 * t36;
	t23 = -t28 * t39 + t29 * t36;
	t21 = 0.1e1 / t23 ^ 2;
	t22 = t28 * t38 + t29 * t35;
	t37 = t22 ^ 2 * t21 + 0.1e1;
	t19 = 0.1e1 / t37;
	t1 = [0, -t28 * t33 * t41, 0, 0, 0; 0, (0.1e1 / t20 * t40 - (-t25 * t29 * t30 * t41 + (t26 - 0.1e1) * t33 * t24) * t33 * t42) / (t30 * t42 + 0.1e1), 0, 0, 0; 0, ((-t28 * t35 + t29 * t38) / t23 - (-t28 * t36 - t29 * t39) * t22 * t21) * t19, t37 * t19, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:36
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (292->20), mult. (182->39), div. (32->9), fcn. (230->13), ass. (0->29)
	t65 = cos(pkin(5));
	t62 = pkin(10) + qJ(2);
	t55 = cos(t62);
	t64 = sin(pkin(5));
	t66 = t55 * t64;
	t52 = atan2(t66, t65);
	t47 = sin(t52);
	t48 = cos(t52);
	t42 = t47 * t66 + t48 * t65;
	t54 = sin(t62);
	t69 = 0.1e1 / t42 ^ 2 * t54 ^ 2;
	t63 = qJ(3) + qJ(4);
	t56 = pkin(5) + t63;
	t57 = pkin(5) - t63;
	t49 = sin(t56) / 0.2e1 - sin(t57) / 0.2e1;
	t59 = cos(t63);
	t46 = -t54 * t49 + t55 * t59;
	t44 = 0.1e1 / t46 ^ 2;
	t50 = cos(t57) / 0.2e1 + cos(t56) / 0.2e1;
	t58 = sin(t63);
	t45 = t54 * t50 + t55 * t58;
	t68 = t45 ^ 2 * t44;
	t60 = t64 ^ 2;
	t51 = 0.1e1 / (0.1e1 + t55 ^ 2 * t60 / t65 ^ 2);
	t67 = t51 / t65;
	t43 = 0.1e1 / t46;
	t39 = 0.1e1 / (0.1e1 + t68);
	t38 = (t46 * t43 + t68) * t39;
	t1 = [0, -t54 * t64 * t67, 0, 0, 0; 0, (0.1e1 / t42 * t66 - (-t48 * t55 * t60 * t67 + (t51 - 0.1e1) * t64 * t47) * t64 * t69) / (t60 * t69 + 0.1e1), 0, 0, 0; 0, ((t55 * t50 - t54 * t58) * t43 - (-t55 * t49 - t54 * t59) * t45 * t44) * t39, t38, t38, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 21:46:36
	% EndTime: 2024-09-27 21:46:37
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (484->20), mult. (224->39), div. (38->9), fcn. (275->13), ass. (0->29)
	t75 = cos(pkin(5));
	t73 = pkin(10) + qJ(2);
	t69 = cos(t73);
	t74 = sin(pkin(5));
	t76 = t69 * t74;
	t62 = atan2(t76, t75);
	t59 = sin(t62);
	t60 = cos(t62);
	t52 = t59 * t76 + t60 * t75;
	t68 = sin(t73);
	t79 = 0.1e1 / t52 ^ 2 * t68 ^ 2;
	t70 = qJ(3) + qJ(4) + qJ(5);
	t64 = pkin(5) + t70;
	t65 = pkin(5) - t70;
	t57 = sin(t64) / 0.2e1 - sin(t65) / 0.2e1;
	t67 = cos(t70);
	t56 = -t68 * t57 + t69 * t67;
	t54 = 0.1e1 / t56 ^ 2;
	t58 = cos(t65) / 0.2e1 + cos(t64) / 0.2e1;
	t66 = sin(t70);
	t55 = t68 * t58 + t69 * t66;
	t78 = t55 ^ 2 * t54;
	t71 = t74 ^ 2;
	t61 = 0.1e1 / (0.1e1 + t69 ^ 2 * t71 / t75 ^ 2);
	t77 = t61 / t75;
	t53 = 0.1e1 / t56;
	t49 = 0.1e1 / (0.1e1 + t78);
	t48 = (t56 * t53 + t78) * t49;
	t1 = [0, -t68 * t74 * t77, 0, 0, 0; 0, (0.1e1 / t52 * t76 - (-t60 * t69 * t71 * t77 + (t61 - 0.1e1) * t74 * t59) * t74 * t79) / (t71 * t79 + 0.1e1), 0, 0, 0; 0, ((t69 * t58 - t68 * t66) * t53 - (-t69 * t57 - t68 * t67) * t55 * t54) * t49, t48, t48, t48;];
	Ja_rot = t1;
end