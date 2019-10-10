% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRP6
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
%   Wie in S6RPPRRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:56
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRP6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (63->18), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->35)
	t40 = sin(qJ(4));
	t41 = sin(qJ(1));
	t43 = cos(qJ(4));
	t49 = t41 * t43;
	t35 = atan2(t49, t40);
	t32 = sin(t35);
	t33 = cos(t35);
	t25 = t32 * t49 + t33 * t40;
	t24 = 0.1e1 / t25 ^ 2;
	t44 = cos(qJ(1));
	t56 = t24 * t44 ^ 2;
	t42 = cos(qJ(5));
	t46 = t44 * t42;
	t39 = sin(qJ(5));
	t51 = t41 * t39;
	t31 = t40 * t46 - t51;
	t29 = 0.1e1 / t31 ^ 2;
	t47 = t44 * t39;
	t50 = t41 * t42;
	t30 = t40 * t47 + t50;
	t55 = t29 * t30;
	t54 = t32 * t40;
	t38 = t43 ^ 2;
	t53 = 0.1e1 / t40 ^ 2 * t38;
	t34 = 0.1e1 / (t41 ^ 2 * t53 + 0.1e1);
	t52 = t41 * t34;
	t48 = t43 * t44;
	t45 = t30 ^ 2 * t29 + 0.1e1;
	t36 = 0.1e1 / t40;
	t28 = 0.1e1 / t31;
	t27 = (-0.1e1 - t53) * t52;
	t26 = 0.1e1 / t45;
	t23 = 0.1e1 / t25;
	t22 = 0.1e1 / (t38 * t56 + 0.1e1);
	t1 = [t36 * t34 * t48, 0, 0, t27, 0, 0; (t23 * t49 + (t33 * t36 * t38 * t52 + (-t34 + 0.1e1) * t43 * t32) * t43 * t56) * t22, 0, 0, (t40 * t23 + (-t41 * t54 + t33 * t43 + (t33 * t49 - t54) * t27) * t43 * t24) * t44 * t22, 0, 0; ((-t40 * t51 + t46) * t28 - (-t40 * t50 - t47) * t55) * t26, 0, 0, (t28 * t39 - t42 * t55) * t26 * t48, t45 * t26, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:18
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (160->25), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->40)
	t54 = cos(qJ(4));
	t68 = t54 ^ 2;
	t51 = sin(qJ(4));
	t53 = cos(qJ(5));
	t55 = cos(qJ(1));
	t57 = t55 * t53;
	t50 = sin(qJ(5));
	t52 = sin(qJ(1));
	t62 = t52 * t50;
	t38 = t51 * t62 - t57;
	t60 = t54 * t50;
	t35 = atan2(-t38, t60);
	t31 = sin(t35);
	t32 = cos(t35);
	t30 = -t31 * t38 + t32 * t60;
	t29 = 0.1e1 / t30 ^ 2;
	t58 = t55 * t50;
	t61 = t52 * t53;
	t41 = t51 * t58 + t61;
	t67 = t29 * t41;
	t65 = t32 * t38;
	t64 = t41 ^ 2 * t29;
	t44 = 0.1e1 / t50;
	t47 = 0.1e1 / t54;
	t63 = t44 * t47;
	t59 = t54 * t55;
	t42 = t51 * t57 - t62;
	t37 = 0.1e1 / t42 ^ 2;
	t56 = t55 ^ 2 * t68 * t37;
	t48 = 0.1e1 / t68;
	t45 = 0.1e1 / t50 ^ 2;
	t40 = t51 * t61 + t58;
	t36 = 0.1e1 / t42;
	t34 = 0.1e1 / (t38 ^ 2 * t48 * t45 + 0.1e1);
	t33 = 0.1e1 / (0.1e1 + t56);
	t28 = 0.1e1 / t30;
	t27 = (-t38 * t44 * t48 * t51 - t52) * t34;
	t26 = 0.1e1 / (0.1e1 + t64);
	t25 = (t38 * t45 * t53 - t40 * t44) * t47 * t34;
	t1 = [-t41 * t34 * t63, 0, 0, t27, t25, 0; (-t38 * t28 - (-t31 + (t63 * t65 + t31) * t34) * t64) * t26, 0, 0, (t27 * t65 * t67 + (t28 * t59 - (-t32 * t51 + (-t27 - t52) * t31 * t54) * t67) * t50) * t26, (t42 * t28 - (t32 * t54 * t53 - t31 * t40 + (-t31 * t60 - t65) * t25) * t67) * t26, 0; (t37 * t40 * t55 - t36 * t52) * t54 * t33, 0, 0, (-t36 * t51 * t55 - t53 * t56) * t33, t41 * t37 * t33 * t59, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end