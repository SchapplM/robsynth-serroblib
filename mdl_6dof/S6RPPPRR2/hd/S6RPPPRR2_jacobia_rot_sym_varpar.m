% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPPRR2
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
%   Wie in S6RPPPRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPPRR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:27:33
	% EndTime: 2019-10-09 23:27:33
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (324->21), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->38)
	t47 = qJ(1) + pkin(9);
	t43 = sin(t47);
	t62 = t43 ^ 2;
	t46 = pkin(10) + qJ(5);
	t42 = sin(t46);
	t44 = cos(t46);
	t45 = cos(t47);
	t53 = t45 * t44;
	t37 = atan2(-t53, t42);
	t35 = sin(t37);
	t36 = cos(t37);
	t28 = -t35 * t53 + t36 * t42;
	t27 = 0.1e1 / t28 ^ 2;
	t61 = t27 * t44;
	t48 = sin(qJ(6));
	t52 = t45 * t48;
	t49 = cos(qJ(6));
	t55 = t43 * t49;
	t34 = t42 * t55 + t52;
	t32 = 0.1e1 / t34 ^ 2;
	t51 = t45 * t49;
	t56 = t43 * t48;
	t33 = t42 * t56 - t51;
	t60 = t32 * t33;
	t59 = t35 * t42;
	t41 = t44 ^ 2;
	t58 = 0.1e1 / t42 ^ 2 * t41;
	t57 = t43 * t44;
	t38 = 0.1e1 / (t45 ^ 2 * t58 + 0.1e1);
	t54 = t45 * t38;
	t50 = t33 ^ 2 * t32 + 0.1e1;
	t39 = 0.1e1 / t42;
	t31 = 0.1e1 / t34;
	t30 = 0.1e1 / t50;
	t29 = (0.1e1 + t58) * t54;
	t26 = 0.1e1 / t28;
	t25 = 0.1e1 / (t62 * t41 * t27 + 0.1e1);
	t1 = [t39 * t38 * t57, 0, 0, 0, t29, 0; (-t26 * t53 + (-t36 * t39 * t41 * t54 + (-t38 + 0.1e1) * t44 * t35) * t62 * t61) * t25, 0, 0, 0, (t42 * t26 + (t45 * t59 + t36 * t44 + (-t36 * t53 - t59) * t29) * t61) * t43 * t25, 0; ((t42 * t52 + t55) * t31 - (t42 * t51 - t56) * t60) * t30, 0, 0, 0, (t31 * t48 - t49 * t60) * t30 * t57, t50 * t30;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end