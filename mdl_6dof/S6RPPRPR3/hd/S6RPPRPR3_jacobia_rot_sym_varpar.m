% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR3
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
%   Wie in S6RPPRPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRPR3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:42
	% EndTime: 2019-10-09 23:37:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:42
	% EndTime: 2019-10-09 23:37:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:42
	% EndTime: 2019-10-09 23:37:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:42
	% EndTime: 2019-10-09 23:37:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:42
	% EndTime: 2019-10-09 23:37:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:42
	% EndTime: 2019-10-09 23:37:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:37:42
	% EndTime: 2019-10-09 23:37:43
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (324->21), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->38)
	t50 = qJ(1) + pkin(9);
	t46 = sin(t50);
	t65 = t46 ^ 2;
	t49 = qJ(4) + pkin(10);
	t45 = sin(t49);
	t47 = cos(t49);
	t48 = cos(t50);
	t56 = t48 * t47;
	t40 = atan2(-t56, t45);
	t38 = sin(t40);
	t39 = cos(t40);
	t31 = -t38 * t56 + t39 * t45;
	t30 = 0.1e1 / t31 ^ 2;
	t64 = t30 * t47;
	t51 = sin(qJ(6));
	t55 = t48 * t51;
	t52 = cos(qJ(6));
	t58 = t46 * t52;
	t37 = t45 * t58 + t55;
	t35 = 0.1e1 / t37 ^ 2;
	t54 = t48 * t52;
	t59 = t46 * t51;
	t36 = t45 * t59 - t54;
	t63 = t35 * t36;
	t62 = t38 * t45;
	t44 = t47 ^ 2;
	t61 = 0.1e1 / t45 ^ 2 * t44;
	t60 = t46 * t47;
	t41 = 0.1e1 / (t48 ^ 2 * t61 + 0.1e1);
	t57 = t48 * t41;
	t53 = t36 ^ 2 * t35 + 0.1e1;
	t42 = 0.1e1 / t45;
	t34 = 0.1e1 / t37;
	t33 = 0.1e1 / t53;
	t32 = (0.1e1 + t61) * t57;
	t29 = 0.1e1 / t31;
	t28 = 0.1e1 / (t65 * t44 * t30 + 0.1e1);
	t1 = [t42 * t41 * t60, 0, 0, t32, 0, 0; (-t29 * t56 + (-t39 * t42 * t44 * t57 + (-t41 + 0.1e1) * t47 * t38) * t65 * t64) * t28, 0, 0, (t45 * t29 + (t48 * t62 + t39 * t47 + (-t39 * t56 - t62) * t32) * t64) * t46 * t28, 0, 0; ((t45 * t55 + t58) * t34 - (t45 * t54 - t59) * t63) * t33, 0, 0, (t34 * t51 - t52 * t63) * t33 * t60, 0, t53 * t33;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end