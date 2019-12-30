% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRPR14
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
%   Wie in S5RPRPR14_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPRPR14_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR14_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:08:08
	% EndTime: 2019-12-29 17:08:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:08:08
	% EndTime: 2019-12-29 17:08:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:08:08
	% EndTime: 2019-12-29 17:08:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:08:08
	% EndTime: 2019-12-29 17:08:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:08:08
	% EndTime: 2019-12-29 17:08:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:08:08
	% EndTime: 2019-12-29 17:08:08
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (215->20), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
	t46 = qJ(3) + pkin(8);
	t44 = sin(t46);
	t45 = cos(t46);
	t50 = cos(qJ(1));
	t54 = t50 * t45;
	t39 = atan2(-t54, t44);
	t37 = sin(t39);
	t38 = cos(t39);
	t30 = -t37 * t54 + t38 * t44;
	t29 = 0.1e1 / t30 ^ 2;
	t48 = sin(qJ(1));
	t62 = t29 * t48 ^ 2;
	t47 = sin(qJ(5));
	t53 = t50 * t47;
	t49 = cos(qJ(5));
	t56 = t48 * t49;
	t36 = t44 * t56 + t53;
	t34 = 0.1e1 / t36 ^ 2;
	t52 = t50 * t49;
	t57 = t48 * t47;
	t35 = t44 * t57 - t52;
	t61 = t34 * t35;
	t60 = t37 * t44;
	t43 = t45 ^ 2;
	t59 = 0.1e1 / t44 ^ 2 * t43;
	t58 = t45 * t48;
	t40 = 0.1e1 / (t50 ^ 2 * t59 + 0.1e1);
	t55 = t50 * t40;
	t51 = t35 ^ 2 * t34 + 0.1e1;
	t41 = 0.1e1 / t44;
	t33 = 0.1e1 / t36;
	t32 = 0.1e1 / t51;
	t31 = (0.1e1 + t59) * t55;
	t28 = 0.1e1 / t30;
	t27 = 0.1e1 / (t43 * t62 + 0.1e1);
	t1 = [t41 * t40 * t58, 0, t31, 0, 0; (-t28 * t54 + (-t38 * t41 * t43 * t55 + (-t40 + 0.1e1) * t45 * t37) * t45 * t62) * t27, 0, (t44 * t28 + (t50 * t60 + t38 * t45 + (-t38 * t54 - t60) * t31) * t45 * t29) * t48 * t27, 0, 0; ((t44 * t53 + t56) * t33 - (t44 * t52 - t57) * t61) * t32, 0, (t33 * t47 - t49 * t61) * t32 * t58, 0, t51 * t32;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end