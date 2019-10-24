% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPPRR2
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
%   Wie in S5PPPRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:17
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PPPRR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPPRR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:40
	% EndTime: 2019-10-24 10:17:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->0), mult. (48->0), div. (5->0), fcn. (63->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 1, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:41
	% EndTime: 2019-10-24 10:17:41
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (186->22), mult. (521->50), div. (35->9), fcn. (739->13), ass. (0->36)
	t55 = sin(pkin(8));
	t61 = sin(qJ(4));
	t68 = t55 * t61;
	t63 = cos(qJ(4));
	t67 = t55 * t63;
	t57 = cos(pkin(9));
	t58 = cos(pkin(8));
	t66 = t57 * t58;
	t54 = sin(pkin(9));
	t59 = cos(pkin(7));
	t65 = t59 * t54;
	t56 = sin(pkin(7));
	t48 = t56 * t54 + t59 * t66;
	t45 = t48 * t63 + t59 * t68;
	t47 = -t56 * t57 + t58 * t65;
	t60 = sin(qJ(5));
	t62 = cos(qJ(5));
	t36 = t45 * t62 + t47 * t60;
	t34 = 0.1e1 / t36 ^ 2;
	t35 = t45 * t60 - t47 * t62;
	t64 = t35 ^ 2 * t34 + 0.1e1;
	t51 = t57 * t67 - t58 * t61;
	t50 = t57 * t68 + t58 * t63;
	t49 = 0.1e1 / t50 ^ 2;
	t46 = t56 * t66 - t65;
	t44 = t48 * t61 - t59 * t67;
	t43 = t46 * t63 + t56 * t68;
	t41 = t46 * t61 - t56 * t67;
	t40 = atan2(-t41, t50);
	t38 = cos(t40);
	t37 = sin(t40);
	t33 = 0.1e1 / t64;
	t32 = -t37 * t41 + t38 * t50;
	t31 = 0.1e1 / t32 ^ 2;
	t29 = (-t43 / t50 + t51 * t41 * t49) / (t41 ^ 2 * t49 + 0.1e1);
	t1 = [0, 0, 0, t29, 0; 0, 0, 0, (t45 / t32 - (-t37 * t43 + t38 * t51 + (-t37 * t50 - t38 * t41) * t29) * t44 * t31) / (t44 ^ 2 * t31 + 0.1e1), 0; 0, 0, 0, (-t60 / t36 + t62 * t35 * t34) * t44 * t33, t64 * t33;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end