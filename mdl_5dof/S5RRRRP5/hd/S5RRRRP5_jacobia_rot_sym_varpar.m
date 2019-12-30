% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRP5
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
%   Wie in S5RRRRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:30
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRRRP5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP5_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP5_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:30
	% EndTime: 2019-12-29 20:30:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:25
	% EndTime: 2019-12-29 20:30:25
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:25
	% EndTime: 2019-12-29 20:30:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:25
	% EndTime: 2019-12-29 20:30:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:25
	% EndTime: 2019-12-29 20:30:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:26
	% EndTime: 2019-12-29 20:30:26
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (638->19), mult. (309->41), div. (101->9), fcn. (474->7), ass. (0->29)
	t60 = cos(qJ(1));
	t61 = t60 ^ 2;
	t55 = qJ(2) + qJ(3) + qJ(4);
	t54 = cos(t55);
	t53 = sin(t55);
	t59 = sin(qJ(1));
	t63 = t59 * t53;
	t47 = atan2(-t63, -t54);
	t45 = sin(t47);
	t46 = cos(t47);
	t43 = -t45 * t63 - t46 * t54;
	t42 = 0.1e1 / t43 ^ 2;
	t68 = t42 * t53;
	t67 = t45 * t54;
	t50 = t53 ^ 2;
	t52 = 0.1e1 / t54 ^ 2;
	t66 = t50 * t52;
	t56 = t59 ^ 2;
	t65 = t56 / t61;
	t48 = 0.1e1 / (t56 * t66 + 0.1e1);
	t64 = t59 * t48;
	t49 = 0.1e1 / (t52 * t65 + 0.1e1);
	t62 = 0.1e1 / t60 * t52 * t49 * t63;
	t51 = 0.1e1 / t54;
	t44 = (0.1e1 + t66) * t64;
	t41 = 0.1e1 / t43;
	t40 = 0.1e1 / (t61 * t50 * t42 + 0.1e1);
	t39 = (t54 * t41 - (-t59 * t67 + t46 * t53 + (-t46 * t63 + t67) * t44) * t68) * t60 * t40;
	t1 = [t60 * t53 * t51 * t48, t44, t44, t44, 0; (-t41 * t63 - (-t46 * t50 * t51 * t64 + (t48 - 0.1e1) * t53 * t45) * t61 * t68) * t40, t39, t39, t39, 0; (-0.1e1 - t65) * t51 * t49, -t62, -t62, -t62, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end