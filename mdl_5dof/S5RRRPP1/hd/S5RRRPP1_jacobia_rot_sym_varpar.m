% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPP1
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
%   Wie in S5RRRPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRRPP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:50:40
	% EndTime: 2019-12-31 20:50:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:50:40
	% EndTime: 2019-12-31 20:50:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:50:40
	% EndTime: 2019-12-31 20:50:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (18->0), mult. (12->0), div. (10->0), fcn. (12->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:50:40
	% EndTime: 2019-12-31 20:50:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:50:40
	% EndTime: 2019-12-31 20:50:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:50:40
	% EndTime: 2019-12-31 20:50:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (384->18), mult. (208->39), div. (72->9), fcn. (333->7), ass. (0->31)
	t51 = qJ(1) + qJ(2);
	t49 = cos(t51);
	t52 = t49 ^ 2;
	t50 = qJ(3) + pkin(8);
	t47 = cos(t50);
	t46 = sin(t50);
	t48 = sin(t51);
	t53 = t48 * t46;
	t38 = atan2(-t53, -t47);
	t36 = sin(t38);
	t37 = cos(t38);
	t32 = -t36 * t53 - t37 * t47;
	t31 = 0.1e1 / t32 ^ 2;
	t58 = t31 * t46;
	t57 = t36 * t47;
	t41 = t46 ^ 2;
	t43 = 0.1e1 / t47 ^ 2;
	t56 = t41 * t43;
	t44 = t48 ^ 2;
	t55 = t44 / t52;
	t39 = 0.1e1 / (t44 * t56 + 0.1e1);
	t54 = t48 * t39;
	t42 = 0.1e1 / t47;
	t40 = 0.1e1 / (t43 * t55 + 0.1e1);
	t35 = t49 * t46 * t42 * t39;
	t34 = (-0.1e1 - t55) * t42 * t40;
	t33 = (0.1e1 + t56) * t54;
	t30 = 0.1e1 / t32;
	t29 = 0.1e1 / (t52 * t41 * t31 + 0.1e1);
	t28 = (-t30 * t53 - (-t37 * t41 * t42 * t54 + (t39 - 0.1e1) * t46 * t36) * t52 * t58) * t29;
	t1 = [t35, t35, t33, 0, 0; t28, t28, (t47 * t30 - (-t48 * t57 + t37 * t46 + (-t37 * t53 + t57) * t33) * t58) * t49 * t29, 0, 0; t34, t34, -0.1e1 / t49 * t43 * t40 * t53, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end