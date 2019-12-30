% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPR8
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
%   Wie in S5RRRPR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRRPR8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:07:57
	% EndTime: 2019-12-29 20:07:58
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (297->18), mult. (231->39), div. (67->9), fcn. (349->7), ass. (0->30)
	t55 = cos(qJ(1));
	t52 = t55 ^ 2;
	t53 = qJ(2) + qJ(3);
	t49 = cos(t53);
	t48 = sin(t53);
	t54 = sin(qJ(1));
	t59 = t54 * t48;
	t42 = atan2(-t59, -t49);
	t40 = sin(t42);
	t41 = cos(t42);
	t38 = -t40 * t59 - t41 * t49;
	t37 = 0.1e1 / t38 ^ 2;
	t65 = t37 * t48;
	t64 = t40 * t49;
	t45 = t48 ^ 2;
	t56 = t49 ^ 2;
	t63 = t45 / t56;
	t62 = t48 * t55;
	t57 = t54 ^ 2;
	t61 = 0.1e1 / t57 * t52;
	t43 = 0.1e1 / (t57 * t63 + 0.1e1);
	t60 = t54 * t43;
	t44 = 0.1e1 / (t56 * t61 + 0.1e1);
	t58 = 0.1e1 / t54 * t44 * t62;
	t46 = 0.1e1 / t49;
	t39 = (0.1e1 + t63) * t60;
	t36 = 0.1e1 / t38;
	t35 = 0.1e1 / (t52 * t45 * t37 + 0.1e1);
	t34 = (t49 * t36 - (-t54 * t64 + t41 * t48 + (-t41 * t59 + t64) * t39) * t65) * t55 * t35;
	t1 = [t46 * t43 * t62, t39, t39, 0, 0; (-t36 * t59 - (-t41 * t45 * t46 * t60 + (t43 - 0.1e1) * t48 * t40) * t52 * t65) * t35, t34, t34, 0, 0; (-0.1e1 - t61) * t49 * t44, -t58, -t58, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:08:03
	% EndTime: 2019-12-29 20:08:03
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (324->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t71 = qJ(2) + qJ(3);
	t69 = sin(t71);
	t70 = cos(t71);
	t73 = sin(qJ(1));
	t81 = t73 * t70;
	t64 = atan2(-t81, t69);
	t62 = sin(t64);
	t63 = cos(t64);
	t55 = -t62 * t81 + t63 * t69;
	t54 = 0.1e1 / t55 ^ 2;
	t75 = cos(qJ(1));
	t87 = t54 * t75 ^ 2;
	t72 = sin(qJ(5));
	t78 = t75 * t72;
	t74 = cos(qJ(5));
	t79 = t73 * t74;
	t61 = t69 * t78 + t79;
	t59 = 0.1e1 / t61 ^ 2;
	t77 = t75 * t74;
	t80 = t73 * t72;
	t60 = -t69 * t77 + t80;
	t86 = t59 * t60;
	t85 = t62 * t69;
	t68 = t70 ^ 2;
	t84 = 0.1e1 / t69 ^ 2 * t68;
	t83 = t70 * t75;
	t65 = 0.1e1 / (t73 ^ 2 * t84 + 0.1e1);
	t82 = t73 * t65;
	t76 = t60 ^ 2 * t59 + 0.1e1;
	t66 = 0.1e1 / t69;
	t58 = 0.1e1 / t61;
	t57 = 0.1e1 / t76;
	t56 = (0.1e1 + t84) * t82;
	t53 = 0.1e1 / t55;
	t52 = 0.1e1 / (t68 * t87 + 0.1e1);
	t51 = (-t58 * t74 - t72 * t86) * t57 * t83;
	t50 = (-t69 * t53 - (t73 * t85 + t63 * t70 + (-t63 * t81 - t85) * t56) * t70 * t54) * t75 * t52;
	t1 = [-t66 * t65 * t83, t56, t56, 0, 0; (-t53 * t81 - (t63 * t66 * t68 * t82 + (t65 - 0.1e1) * t70 * t62) * t70 * t87) * t52, t50, t50, 0, 0; ((t69 * t79 + t78) * t58 - (-t69 * t80 + t77) * t86) * t57, t51, t51, 0, t76 * t57;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end