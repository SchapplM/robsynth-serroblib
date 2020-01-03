% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPP3
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
%   Wie in S5RRRPP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRRPP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_jacobia_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:54:02
	% EndTime: 2019-12-31 20:54:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:54:02
	% EndTime: 2019-12-31 20:54:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:54:02
	% EndTime: 2019-12-31 20:54:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (18->0), mult. (12->0), div. (10->0), fcn. (12->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:54:02
	% EndTime: 2019-12-31 20:54:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:54:02
	% EndTime: 2019-12-31 20:54:02
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (223->17), mult. (216->38), div. (62->9), fcn. (332->7), ass. (0->31)
	t46 = qJ(1) + qJ(2);
	t42 = cos(t46);
	t40 = t42 ^ 2;
	t48 = cos(qJ(3));
	t41 = sin(t46);
	t47 = sin(qJ(3));
	t53 = t41 * t47;
	t37 = atan2(-t53, -t48);
	t34 = sin(t37);
	t35 = cos(t37);
	t31 = -t34 * t53 - t35 * t48;
	t30 = 0.1e1 / t31 ^ 2;
	t57 = t30 * t47;
	t56 = t34 * t48;
	t49 = t41 ^ 2;
	t55 = 0.1e1 / t49 * t40;
	t43 = t47 ^ 2;
	t50 = t48 ^ 2;
	t51 = t43 / t50;
	t38 = 0.1e1 / (t49 * t51 + 0.1e1);
	t54 = t41 * t38;
	t52 = t42 * t47;
	t44 = 0.1e1 / t48;
	t36 = 0.1e1 / (t50 * t55 + 0.1e1);
	t33 = t44 * t38 * t52;
	t32 = (0.1e1 + t51) * t54;
	t29 = 0.1e1 / t31;
	t28 = (-0.1e1 - t55) * t48 * t36;
	t27 = 0.1e1 / (t40 * t43 * t30 + 0.1e1);
	t26 = (-t29 * t53 - (-t35 * t43 * t44 * t54 + (t38 - 0.1e1) * t47 * t34) * t40 * t57) * t27;
	t1 = [t33, t33, t32, 0, 0; t26, t26, (t48 * t29 - (-t41 * t56 + t35 * t47 + (-t35 * t53 + t56) * t32) * t57) * t42 * t27, 0, 0; t28, t28, -0.1e1 / t41 * t36 * t52, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:54:02
	% EndTime: 2019-12-31 20:54:02
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (191->18), mult. (216->38), div. (62->9), fcn. (332->7), ass. (0->31)
	t46 = qJ(1) + qJ(2);
	t42 = cos(t46);
	t40 = t42 ^ 2;
	t47 = sin(qJ(3));
	t41 = sin(t46);
	t48 = cos(qJ(3));
	t54 = t41 * t48;
	t37 = atan2(-t54, t47);
	t35 = sin(t37);
	t36 = cos(t37);
	t32 = -t35 * t54 + t36 * t47;
	t31 = 0.1e1 / t32 ^ 2;
	t58 = t31 * t48;
	t57 = t35 * t47;
	t49 = t41 ^ 2;
	t56 = 0.1e1 / t49 * t40;
	t45 = t48 ^ 2;
	t50 = t47 ^ 2;
	t52 = 0.1e1 / t50 * t45;
	t38 = 0.1e1 / (t49 * t52 + 0.1e1);
	t55 = t41 * t38;
	t53 = t42 * t48;
	t43 = 0.1e1 / t47;
	t51 = t43 * t38 * t53;
	t34 = 0.1e1 / (t50 * t56 + 0.1e1);
	t33 = (0.1e1 + t52) * t55;
	t30 = 0.1e1 / t32;
	t29 = (0.1e1 + t56) * t47 * t34;
	t28 = 0.1e1 / (t40 * t45 * t31 + 0.1e1);
	t27 = (-t30 * t54 - (t36 * t43 * t45 * t55 + (t38 - 0.1e1) * t48 * t35) * t40 * t58) * t28;
	t1 = [-t51, -t51, t33, 0, 0; t27, t27, (-t47 * t30 - (t41 * t57 + t36 * t48 + (-t36 * t54 - t57) * t33) * t58) * t42 * t28, 0, 0; t29, t29, -0.1e1 / t41 * t34 * t53, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end