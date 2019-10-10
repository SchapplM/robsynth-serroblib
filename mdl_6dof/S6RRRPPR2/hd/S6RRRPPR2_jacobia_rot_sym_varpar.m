% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR2
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
%   Wie in S6RRRPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (467->18), mult. (231->39), div. (67->9), fcn. (349->7), ass. (0->30)
	t61 = cos(qJ(1));
	t59 = t61 ^ 2;
	t56 = qJ(2) + qJ(3) + pkin(10);
	t55 = cos(t56);
	t54 = sin(t56);
	t60 = sin(qJ(1));
	t65 = t60 * t54;
	t48 = atan2(-t65, -t55);
	t46 = sin(t48);
	t47 = cos(t48);
	t44 = -t46 * t65 - t47 * t55;
	t43 = 0.1e1 / t44 ^ 2;
	t71 = t43 * t54;
	t70 = t46 * t55;
	t51 = t54 ^ 2;
	t62 = t55 ^ 2;
	t69 = t51 / t62;
	t68 = t54 * t61;
	t63 = t60 ^ 2;
	t67 = 0.1e1 / t63 * t59;
	t49 = 0.1e1 / (t63 * t69 + 0.1e1);
	t66 = t60 * t49;
	t50 = 0.1e1 / (t62 * t67 + 0.1e1);
	t64 = 0.1e1 / t60 * t50 * t68;
	t52 = 0.1e1 / t55;
	t45 = (0.1e1 + t69) * t66;
	t42 = 0.1e1 / t44;
	t41 = 0.1e1 / (t59 * t51 * t43 + 0.1e1);
	t40 = (t55 * t42 - (-t60 * t70 + t47 * t54 + (-t47 * t65 + t70) * t45) * t71) * t61 * t41;
	t1 = [t52 * t49 * t68, t45, t45, 0, 0, 0; (-t42 * t65 - (-t47 * t51 * t52 * t66 + (t49 - 0.1e1) * t54 * t46) * t59 * t71) * t41, t40, t40, 0, 0, 0; (-0.1e1 - t67) * t55 * t50, -t64, -t64, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (520->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t77 = qJ(2) + qJ(3) + pkin(10);
	t75 = sin(t77);
	t76 = cos(t77);
	t79 = sin(qJ(1));
	t87 = t79 * t76;
	t70 = atan2(-t87, t75);
	t64 = sin(t70);
	t65 = cos(t70);
	t61 = -t64 * t87 + t65 * t75;
	t60 = 0.1e1 / t61 ^ 2;
	t81 = cos(qJ(1));
	t93 = t60 * t81 ^ 2;
	t92 = t64 * t75;
	t78 = sin(qJ(6));
	t84 = t81 * t78;
	t80 = cos(qJ(6));
	t85 = t79 * t80;
	t69 = t75 * t84 + t85;
	t67 = 0.1e1 / t69 ^ 2;
	t83 = t81 * t80;
	t86 = t79 * t78;
	t68 = -t75 * t83 + t86;
	t91 = t67 * t68;
	t74 = t76 ^ 2;
	t90 = 0.1e1 / t75 ^ 2 * t74;
	t89 = t76 * t81;
	t71 = 0.1e1 / (t79 ^ 2 * t90 + 0.1e1);
	t88 = t79 * t71;
	t82 = t68 ^ 2 * t67 + 0.1e1;
	t72 = 0.1e1 / t75;
	t66 = 0.1e1 / t69;
	t63 = 0.1e1 / t82;
	t62 = (0.1e1 + t90) * t88;
	t59 = 0.1e1 / t61;
	t58 = 0.1e1 / (t74 * t93 + 0.1e1);
	t57 = (-t66 * t80 - t78 * t91) * t63 * t89;
	t56 = (-t75 * t59 - (t79 * t92 + t65 * t76 + (-t65 * t87 - t92) * t62) * t76 * t60) * t81 * t58;
	t1 = [-t72 * t71 * t89, t62, t62, 0, 0, 0; (-t59 * t87 - (t65 * t72 * t74 * t88 + (t71 - 0.1e1) * t76 * t64) * t76 * t93) * t58, t56, t56, 0, 0, 0; ((t75 * t85 + t84) * t66 - (-t75 * t86 + t83) * t91) * t63, t57, t57, 0, 0, t82 * t63;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end