% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRPR2
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
%   Wie in S5PRRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:30
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRRPR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:16
	% EndTime: 2019-10-24 10:30:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:16
	% EndTime: 2019-10-24 10:30:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:15
	% EndTime: 2019-10-24 10:30:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:16
	% EndTime: 2019-10-24 10:30:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (30->0), mult. (12->0), div. (10->0), fcn. (12->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 1, 1, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:15
	% EndTime: 2019-10-24 10:30:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:16
	% EndTime: 2019-10-24 10:30:16
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (310->15), mult. (205->35), div. (45->9), fcn. (315->9), ass. (0->27)
	t52 = cos(pkin(9));
	t48 = pkin(8) + qJ(2) + qJ(3);
	t46 = sin(t48);
	t51 = sin(pkin(9));
	t58 = t46 * t51;
	t44 = atan2(-t58, -t52);
	t42 = sin(t44);
	t43 = cos(t44);
	t36 = -t42 * t58 - t43 * t52;
	t47 = cos(t48);
	t60 = 0.1e1 / t36 ^ 2 * t47 ^ 2;
	t49 = t51 ^ 2;
	t45 = 0.1e1 / (0.1e1 + t46 ^ 2 * t49 / t52 ^ 2);
	t59 = t45 / t52;
	t53 = sin(qJ(5));
	t57 = t52 * t53;
	t54 = cos(qJ(5));
	t56 = t52 * t54;
	t41 = t46 * t53 + t47 * t56;
	t39 = 0.1e1 / t41 ^ 2;
	t40 = -t46 * t54 + t47 * t57;
	t55 = t40 ^ 2 * t39 + 0.1e1;
	t38 = t47 * t51 * t59;
	t37 = 0.1e1 / t55;
	t33 = ((-t46 * t57 - t47 * t54) / t41 - (-t46 * t56 + t47 * t53) * t40 * t39) * t37;
	t32 = (-0.1e1 / t36 * t58 - (-t43 * t46 * t49 * t59 + (t45 - 0.1e1) * t51 * t42) * t51 * t60) / (t49 * t60 + 0.1e1);
	t1 = [0, t38, t38, 0, 0; 0, t32, t32, 0, 0; 0, t33, t33, 0, t55 * t37;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end