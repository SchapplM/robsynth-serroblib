% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RPRRPP7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:20
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:20:06
	% EndTime: 2019-10-10 01:20:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:20:06
	% EndTime: 2019-10-10 01:20:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:20:06
	% EndTime: 2019-10-10 01:20:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:20:06
	% EndTime: 2019-10-10 01:20:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:20:07
	% EndTime: 2019-10-10 01:20:07
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (86->19), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->36)
	t40 = sin(qJ(1));
	t56 = t40 ^ 2;
	t39 = sin(qJ(3));
	t42 = cos(qJ(3));
	t43 = cos(qJ(1));
	t45 = t43 * t42;
	t33 = atan2(-t45, t39);
	t31 = sin(t33);
	t32 = cos(t33);
	t24 = -t31 * t45 + t32 * t39;
	t23 = 0.1e1 / t24 ^ 2;
	t55 = t23 * t42;
	t38 = sin(qJ(4));
	t47 = t43 * t38;
	t41 = cos(qJ(4));
	t50 = t40 * t41;
	t30 = t39 * t50 + t47;
	t28 = 0.1e1 / t30 ^ 2;
	t46 = t43 * t41;
	t51 = t40 * t38;
	t29 = t39 * t51 - t46;
	t54 = t28 * t29;
	t53 = t31 * t39;
	t37 = t42 ^ 2;
	t52 = 0.1e1 / t39 ^ 2 * t37;
	t49 = t40 * t42;
	t34 = 0.1e1 / (t43 ^ 2 * t52 + 0.1e1);
	t48 = t43 * t34;
	t44 = t29 ^ 2 * t28 + 0.1e1;
	t35 = 0.1e1 / t39;
	t27 = 0.1e1 / t30;
	t26 = (0.1e1 + t52) * t48;
	t25 = 0.1e1 / t44;
	t22 = 0.1e1 / t24;
	t21 = 0.1e1 / (t56 * t37 * t23 + 0.1e1);
	t1 = [t35 * t34 * t49, 0, t26, 0, 0, 0; (-t22 * t45 + (-t32 * t35 * t37 * t48 + (-t34 + 0.1e1) * t42 * t31) * t56 * t55) * t21, 0, (t39 * t22 + (t43 * t53 + t32 * t42 + (-t32 * t45 - t53) * t26) * t55) * t40 * t21, 0, 0, 0; ((t39 * t47 + t50) * t27 - (t39 * t46 - t51) * t54) * t25, 0, (t27 * t38 - t41 * t54) * t25 * t49, t44 * t25, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:20:07
	% EndTime: 2019-10-10 01:20:07
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (160->24), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->40)
	t54 = cos(qJ(3));
	t68 = t54 ^ 2;
	t51 = sin(qJ(3));
	t50 = sin(qJ(4));
	t55 = cos(qJ(1));
	t58 = t55 * t50;
	t52 = sin(qJ(1));
	t53 = cos(qJ(4));
	t61 = t52 * t53;
	t41 = t51 * t58 + t61;
	t59 = t54 * t50;
	t36 = atan2(t41, t59);
	t32 = sin(t36);
	t33 = cos(t36);
	t31 = t32 * t41 + t33 * t59;
	t30 = 0.1e1 / t31 ^ 2;
	t57 = t55 * t53;
	t62 = t52 * t50;
	t39 = t51 * t62 - t57;
	t67 = t30 * t39;
	t65 = t33 * t41;
	t64 = t39 ^ 2 * t30;
	t44 = 0.1e1 / t50;
	t48 = 0.1e1 / t54;
	t63 = t44 * t48;
	t60 = t52 * t54;
	t40 = t51 * t61 + t58;
	t38 = 0.1e1 / t40 ^ 2;
	t56 = t52 ^ 2 * t68 * t38;
	t49 = 0.1e1 / t68;
	t45 = 0.1e1 / t50 ^ 2;
	t42 = t51 * t57 - t62;
	t37 = 0.1e1 / t40;
	t35 = 0.1e1 / (t41 ^ 2 * t49 * t45 + 0.1e1);
	t34 = 0.1e1 / (0.1e1 + t56);
	t29 = 0.1e1 / t31;
	t28 = (t41 * t44 * t49 * t51 + t55) * t35;
	t27 = 0.1e1 / (0.1e1 + t64);
	t26 = (-t41 * t45 * t53 + t42 * t44) * t48 * t35;
	t1 = [-t39 * t35 * t63, 0, t28, t26, 0, 0; (t41 * t29 - (-t32 + (-t63 * t65 + t32) * t35) * t64) * t27, 0, (-t28 * t65 * t67 + (t29 * t60 - (-t33 * t51 + (-t28 + t55) * t54 * t32) * t67) * t50) * t27, (t40 * t29 - (t33 * t54 * t53 + t32 * t42 + (-t32 * t59 + t65) * t26) * t67) * t27, 0, 0; (-t38 * t42 * t52 + t37 * t55) * t54 * t34, 0, (-t37 * t51 * t52 - t53 * t56) * t34, t39 * t38 * t34 * t60, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:20:07
	% EndTime: 2019-10-10 01:20:07
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (87->20), mult. (227->55), div. (53->9), fcn. (337->9), ass. (0->36)
	t46 = sin(qJ(1));
	t62 = t46 ^ 2;
	t45 = sin(qJ(3));
	t48 = cos(qJ(3));
	t49 = cos(qJ(1));
	t50 = t49 * t48;
	t39 = atan2(t50, -t45);
	t37 = sin(t39);
	t38 = cos(t39);
	t29 = t37 * t50 - t38 * t45;
	t28 = 0.1e1 / t29 ^ 2;
	t61 = t28 * t48;
	t44 = sin(qJ(4));
	t52 = t49 * t44;
	t47 = cos(qJ(4));
	t55 = t46 * t47;
	t36 = t45 * t55 + t52;
	t34 = 0.1e1 / t36 ^ 2;
	t51 = t49 * t47;
	t56 = t46 * t44;
	t35 = -t45 * t56 + t51;
	t60 = t35 ^ 2 * t34;
	t59 = t34 * t35;
	t58 = t37 * t45;
	t43 = t48 ^ 2;
	t57 = 0.1e1 / t45 ^ 2 * t43;
	t54 = t46 * t48;
	t40 = 0.1e1 / (t49 ^ 2 * t57 + 0.1e1);
	t53 = t49 * t40;
	t41 = 0.1e1 / t45;
	t33 = 0.1e1 / t36;
	t31 = (0.1e1 + t57) * t53;
	t30 = 0.1e1 / (0.1e1 + t60);
	t27 = 0.1e1 / t29;
	t26 = 0.1e1 / (t62 * t43 * t28 + 0.1e1);
	t1 = [t41 * t40 * t54, 0, t31, 0, 0, 0; (t27 * t50 - (t38 * t41 * t43 * t53 + (t40 - 0.1e1) * t48 * t37) * t62 * t61) * t26, 0, (-t45 * t27 - (-t49 * t58 - t38 * t48 + (t38 * t50 + t58) * t31) * t61) * t46 * t26, 0, 0, 0; ((-t45 * t52 - t55) * t33 - (t45 * t51 - t56) * t59) * t30, 0, (-t33 * t44 - t47 * t59) * t30 * t54, (-t33 * t36 - t60) * t30, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end