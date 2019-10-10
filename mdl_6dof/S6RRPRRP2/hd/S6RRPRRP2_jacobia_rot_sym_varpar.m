% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP2
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
%   Wie in S6RRPRRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:31
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (554->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t71 = qJ(2) + pkin(10) + qJ(4);
	t70 = cos(t71);
	t69 = sin(t71);
	t73 = sin(qJ(1));
	t81 = t73 * t69;
	t62 = atan2(-t81, -t70);
	t58 = sin(t62);
	t59 = cos(t62);
	t55 = -t58 * t81 - t59 * t70;
	t54 = 0.1e1 / t55 ^ 2;
	t75 = cos(qJ(1));
	t87 = t54 * t75 ^ 2;
	t86 = t58 * t70;
	t74 = cos(qJ(5));
	t77 = t75 * t74;
	t72 = sin(qJ(5));
	t80 = t73 * t72;
	t64 = t70 * t77 + t80;
	t61 = 0.1e1 / t64 ^ 2;
	t78 = t75 * t72;
	t79 = t73 * t74;
	t63 = t70 * t78 - t79;
	t85 = t61 * t63;
	t66 = t69 ^ 2;
	t84 = t66 / t70 ^ 2;
	t83 = t69 * t75;
	t65 = 0.1e1 / (t73 ^ 2 * t84 + 0.1e1);
	t82 = t73 * t65;
	t76 = t63 ^ 2 * t61 + 0.1e1;
	t67 = 0.1e1 / t70;
	t60 = 0.1e1 / t64;
	t57 = 0.1e1 / t76;
	t56 = (0.1e1 + t84) * t82;
	t53 = 0.1e1 / t55;
	t52 = 0.1e1 / (t66 * t87 + 0.1e1);
	t51 = (-t60 * t72 + t74 * t85) * t57 * t83;
	t50 = (t70 * t53 - (-t73 * t86 + t59 * t69 + (-t59 * t81 + t86) * t56) * t69 * t54) * t75 * t52;
	t1 = [t67 * t65 * t83, t56, 0, t56, 0, 0; (-t53 * t81 - (-t59 * t66 * t67 * t82 + (t65 - 0.1e1) * t69 * t58) * t69 * t87) * t52, t50, 0, t50, 0, 0; ((-t70 * t80 - t77) * t60 - (-t70 * t79 + t78) * t85) * t57, t51, 0, t51, t76 * t57, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:31:51
	% EndTime: 2019-10-10 10:31:51
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (740->27), mult. (688->69), div. (136->11), fcn. (1034->9), ass. (0->43)
	t82 = qJ(2) + pkin(10) + qJ(4);
	t78 = sin(t82);
	t102 = t78 ^ 2;
	t79 = cos(t82);
	t88 = cos(qJ(5));
	t89 = cos(qJ(1));
	t91 = t89 * t88;
	t86 = sin(qJ(5));
	t87 = sin(qJ(1));
	t94 = t87 * t86;
	t70 = t79 * t94 + t91;
	t96 = t78 * t86;
	t67 = atan2(-t70, t96);
	t63 = sin(t67);
	t64 = cos(t67);
	t62 = -t63 * t70 + t64 * t96;
	t61 = 0.1e1 / t62 ^ 2;
	t92 = t89 * t86;
	t93 = t87 * t88;
	t73 = t79 * t92 - t93;
	t101 = t61 * t73;
	t99 = t64 * t70;
	t98 = t73 ^ 2 * t61;
	t76 = 0.1e1 / t78;
	t83 = 0.1e1 / t86;
	t97 = t76 * t83;
	t95 = t78 * t89;
	t74 = t79 * t91 + t94;
	t69 = 0.1e1 / t74 ^ 2;
	t90 = t89 ^ 2 * t102 * t69;
	t84 = 0.1e1 / t86 ^ 2;
	t77 = 0.1e1 / t102;
	t72 = t79 * t93 - t92;
	t68 = 0.1e1 / t74;
	t66 = 0.1e1 / (t70 ^ 2 * t77 * t84 + 0.1e1);
	t65 = 0.1e1 / (0.1e1 + t90);
	t60 = 0.1e1 / t62;
	t59 = (t70 * t77 * t79 * t83 + t87) * t66;
	t58 = 0.1e1 / (0.1e1 + t98);
	t57 = (t70 * t84 * t88 - t72 * t83) * t76 * t66;
	t56 = (-t68 * t79 * t89 - t88 * t90) * t65;
	t55 = (t59 * t99 * t101 + (-t60 * t95 - (t64 * t79 + (-t59 + t87) * t78 * t63) * t101) * t86) * t58;
	t1 = [-t73 * t66 * t97, t59, 0, t59, t57, 0; (-t70 * t60 - (-t63 + (t97 * t99 + t63) * t66) * t98) * t58, t55, 0, t55, (t74 * t60 - (t64 * t78 * t88 - t63 * t72 + (-t63 * t96 - t99) * t57) * t101) * t58, 0; (-t69 * t72 * t89 + t68 * t87) * t78 * t65, t56, 0, t56, -t73 * t69 * t65 * t95, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end