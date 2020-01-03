% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRR6
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
%   Wie in S5RPRRR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RPRRR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:02:19
	% EndTime: 2019-12-31 19:02:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:02:19
	% EndTime: 2019-12-31 19:02:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:02:19
	% EndTime: 2019-12-31 19:02:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:02:19
	% EndTime: 2019-12-31 19:02:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:02:19
	% EndTime: 2019-12-31 19:02:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:02:20
	% EndTime: 2019-12-31 19:02:20
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (517->22), mult. (332->57), div. (79->9), fcn. (489->9), ass. (0->38)
	t67 = qJ(1) + pkin(9);
	t64 = cos(t67);
	t81 = t64 ^ 2;
	t68 = qJ(3) + qJ(4);
	t66 = cos(t68);
	t63 = sin(t67);
	t65 = sin(t68);
	t75 = t63 * t65;
	t58 = atan2(-t75, -t66);
	t56 = sin(t58);
	t57 = cos(t58);
	t49 = -t56 * t75 - t57 * t66;
	t48 = 0.1e1 / t49 ^ 2;
	t80 = t48 * t65;
	t69 = sin(qJ(5));
	t70 = cos(qJ(5));
	t72 = t66 * t70;
	t55 = t63 * t69 + t64 * t72;
	t53 = 0.1e1 / t55 ^ 2;
	t73 = t66 * t69;
	t54 = -t63 * t70 + t64 * t73;
	t79 = t53 * t54;
	t78 = t56 * t66;
	t60 = t65 ^ 2;
	t77 = t60 / t66 ^ 2;
	t59 = 0.1e1 / (t63 ^ 2 * t77 + 0.1e1);
	t76 = t63 * t59;
	t74 = t64 * t65;
	t71 = t54 ^ 2 * t53 + 0.1e1;
	t61 = 0.1e1 / t66;
	t52 = 0.1e1 / t55;
	t51 = 0.1e1 / t71;
	t50 = (0.1e1 + t77) * t76;
	t47 = 0.1e1 / t49;
	t46 = 0.1e1 / (t81 * t60 * t48 + 0.1e1);
	t45 = (-t52 * t69 + t70 * t79) * t51 * t74;
	t44 = (t66 * t47 - (-t63 * t78 + t57 * t65 + (-t57 * t75 + t78) * t50) * t80) * t64 * t46;
	t1 = [t61 * t59 * t74, 0, t50, t50, 0; (-t47 * t75 - (-t57 * t60 * t61 * t76 + (t59 - 0.1e1) * t65 * t56) * t81 * t80) * t46, 0, t44, t44, 0; ((-t63 * t73 - t64 * t70) * t52 - (-t63 * t72 + t64 * t69) * t79) * t51, 0, t45, t45, t71 * t51;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end