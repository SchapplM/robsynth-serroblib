% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPP4
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
%   Wie in S5RRRPP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRRPP4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:42:50
	% EndTime: 2019-12-29 19:42:50
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:42:55
	% EndTime: 2019-12-29 19:42:55
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:42:55
	% EndTime: 2019-12-29 19:42:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:42:55
	% EndTime: 2019-12-29 19:42:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:42:55
	% EndTime: 2019-12-29 19:42:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:42:47
	% EndTime: 2019-12-29 19:42:47
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (471->18), mult. (227->41), div. (75->9), fcn. (351->7), ass. (0->29)
	t59 = cos(qJ(1));
	t60 = t59 ^ 2;
	t54 = qJ(2) + qJ(3) + pkin(8);
	t53 = cos(t54);
	t52 = sin(t54);
	t58 = sin(qJ(1));
	t62 = t58 * t52;
	t46 = atan2(-t62, -t53);
	t44 = sin(t46);
	t45 = cos(t46);
	t42 = -t44 * t62 - t45 * t53;
	t41 = 0.1e1 / t42 ^ 2;
	t67 = t41 * t52;
	t66 = t44 * t53;
	t49 = t52 ^ 2;
	t51 = 0.1e1 / t53 ^ 2;
	t65 = t49 * t51;
	t55 = t58 ^ 2;
	t64 = t55 / t60;
	t47 = 0.1e1 / (t55 * t65 + 0.1e1);
	t63 = t58 * t47;
	t48 = 0.1e1 / (t51 * t64 + 0.1e1);
	t61 = 0.1e1 / t59 * t51 * t48 * t62;
	t50 = 0.1e1 / t53;
	t43 = (0.1e1 + t65) * t63;
	t40 = 0.1e1 / t42;
	t39 = 0.1e1 / (t60 * t49 * t41 + 0.1e1);
	t38 = (t53 * t40 - (-t58 * t66 + t45 * t52 + (-t45 * t62 + t66) * t43) * t67) * t59 * t39;
	t1 = [t59 * t52 * t50 * t47, t43, t43, 0, 0; (-t40 * t62 - (-t45 * t49 * t50 * t63 + (t47 - 0.1e1) * t52 * t44) * t60 * t67) * t39, t38, t38, 0, 0; (-0.1e1 - t64) * t50 * t48, -t61, -t61, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end