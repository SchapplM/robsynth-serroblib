% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPPR2
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
%   Wie in S5RRPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRPPR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
Ja_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:21
	% EndTime: 2022-01-20 10:06:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:21
	% EndTime: 2022-01-20 10:06:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:21
	% EndTime: 2022-01-20 10:06:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (18->0), mult. (12->0), div. (10->0), fcn. (12->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:21
	% EndTime: 2022-01-20 10:06:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->0), mult. (12->0), div. (10->0), fcn. (12->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:21
	% EndTime: 2022-01-20 10:06:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 10:06:21
	% EndTime: 2022-01-20 10:06:21
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (310->15), mult. (205->35), div. (45->9), fcn. (315->9), ass. (0->27)
	t54 = cos(pkin(9));
	t50 = qJ(1) + qJ(2) + pkin(8);
	t48 = sin(t50);
	t53 = sin(pkin(9));
	t60 = t48 * t53;
	t46 = atan2(-t60, -t54);
	t44 = sin(t46);
	t45 = cos(t46);
	t38 = -t44 * t60 - t45 * t54;
	t49 = cos(t50);
	t62 = 0.1e1 / t38 ^ 2 * t49 ^ 2;
	t51 = t53 ^ 2;
	t47 = 0.1e1 / (0.1e1 + t48 ^ 2 * t51 / t54 ^ 2);
	t61 = t47 / t54;
	t55 = sin(qJ(5));
	t59 = t54 * t55;
	t56 = cos(qJ(5));
	t58 = t54 * t56;
	t43 = t48 * t55 + t49 * t58;
	t41 = 0.1e1 / t43 ^ 2;
	t42 = -t48 * t56 + t49 * t59;
	t57 = t42 ^ 2 * t41 + 0.1e1;
	t40 = t49 * t53 * t61;
	t39 = 0.1e1 / t57;
	t35 = ((-t48 * t59 - t49 * t56) / t43 - (-t48 * t58 + t49 * t55) * t42 * t41) * t39;
	t34 = (-0.1e1 / t38 * t60 - (-t45 * t48 * t51 * t61 + (t47 - 0.1e1) * t53 * t44) * t53 * t62) / (t51 * t62 + 0.1e1);
	t1 = [t40, t40, 0, 0, 0; t34, t34, 0, 0, 0; t35, t35, 0, 0, t57 * t39;];
	Ja_rot = t1;
end