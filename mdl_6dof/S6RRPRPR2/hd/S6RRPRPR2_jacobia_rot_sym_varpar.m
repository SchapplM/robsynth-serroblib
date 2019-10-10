% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR2
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
%   Wie in S6RRPRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:06
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:08
	% EndTime: 2019-10-10 10:06:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:08
	% EndTime: 2019-10-10 10:06:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:08
	% EndTime: 2019-10-10 10:06:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:08
	% EndTime: 2019-10-10 10:06:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:08
	% EndTime: 2019-10-10 10:06:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:08
	% EndTime: 2019-10-10 10:06:09
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (467->18), mult. (231->39), div. (67->9), fcn. (349->7), ass. (0->30)
	t58 = cos(qJ(1));
	t56 = t58 ^ 2;
	t53 = qJ(2) + pkin(10) + qJ(4);
	t52 = cos(t53);
	t51 = sin(t53);
	t57 = sin(qJ(1));
	t62 = t57 * t51;
	t45 = atan2(-t62, -t52);
	t43 = sin(t45);
	t44 = cos(t45);
	t41 = -t43 * t62 - t44 * t52;
	t40 = 0.1e1 / t41 ^ 2;
	t68 = t40 * t51;
	t67 = t43 * t52;
	t48 = t51 ^ 2;
	t59 = t52 ^ 2;
	t66 = t48 / t59;
	t65 = t51 * t58;
	t60 = t57 ^ 2;
	t64 = 0.1e1 / t60 * t56;
	t46 = 0.1e1 / (t60 * t66 + 0.1e1);
	t63 = t57 * t46;
	t47 = 0.1e1 / (t59 * t64 + 0.1e1);
	t61 = 0.1e1 / t57 * t47 * t65;
	t49 = 0.1e1 / t52;
	t42 = (0.1e1 + t66) * t63;
	t39 = 0.1e1 / t41;
	t38 = 0.1e1 / (t56 * t48 * t40 + 0.1e1);
	t37 = (t52 * t39 - (-t57 * t67 + t44 * t51 + (-t44 * t62 + t67) * t42) * t68) * t58 * t38;
	t1 = [t49 * t46 * t65, t42, 0, t42, 0, 0; (-t39 * t62 - (-t44 * t48 * t49 * t63 + (t46 - 0.1e1) * t51 * t43) * t56 * t68) * t38, t37, 0, t37, 0, 0; (-0.1e1 - t64) * t52 * t47, -t61, 0, -t61, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:09
	% EndTime: 2019-10-10 10:06:09
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (520->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t73 = qJ(2) + pkin(10) + qJ(4);
	t71 = sin(t73);
	t72 = cos(t73);
	t75 = sin(qJ(1));
	t83 = t75 * t72;
	t66 = atan2(-t83, t71);
	t60 = sin(t66);
	t61 = cos(t66);
	t57 = -t60 * t83 + t61 * t71;
	t56 = 0.1e1 / t57 ^ 2;
	t77 = cos(qJ(1));
	t89 = t56 * t77 ^ 2;
	t88 = t60 * t71;
	t74 = sin(qJ(6));
	t80 = t77 * t74;
	t76 = cos(qJ(6));
	t81 = t75 * t76;
	t65 = t71 * t80 + t81;
	t63 = 0.1e1 / t65 ^ 2;
	t79 = t77 * t76;
	t82 = t75 * t74;
	t64 = -t71 * t79 + t82;
	t87 = t63 * t64;
	t70 = t72 ^ 2;
	t86 = 0.1e1 / t71 ^ 2 * t70;
	t85 = t72 * t77;
	t67 = 0.1e1 / (t75 ^ 2 * t86 + 0.1e1);
	t84 = t75 * t67;
	t78 = t64 ^ 2 * t63 + 0.1e1;
	t68 = 0.1e1 / t71;
	t62 = 0.1e1 / t65;
	t59 = 0.1e1 / t78;
	t58 = (0.1e1 + t86) * t84;
	t55 = 0.1e1 / t57;
	t54 = 0.1e1 / (t70 * t89 + 0.1e1);
	t53 = (-t62 * t76 - t74 * t87) * t59 * t85;
	t52 = (-t71 * t55 - (t75 * t88 + t61 * t72 + (-t61 * t83 - t88) * t58) * t72 * t56) * t77 * t54;
	t1 = [-t68 * t67 * t85, t58, 0, t58, 0, 0; (-t55 * t83 - (t61 * t68 * t70 * t84 + (t67 - 0.1e1) * t72 * t60) * t72 * t89) * t54, t52, 0, t52, 0, 0; ((t71 * t81 + t80) * t62 - (-t71 * t82 + t79) * t87) * t59, t53, 0, t53, 0, t78 * t59;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end