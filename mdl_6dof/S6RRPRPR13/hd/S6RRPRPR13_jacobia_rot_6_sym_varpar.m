% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR13
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR13_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:50
% EndTime: 2019-02-26 21:44:50
% DurationCPUTime: 0.15s
% Computational Cost: add. (525->37), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->55)
t78 = cos(pkin(6));
t83 = cos(qJ(2));
t84 = cos(qJ(1));
t89 = t84 * t83;
t80 = sin(qJ(2));
t81 = sin(qJ(1));
t92 = t81 * t80;
t70 = -t78 * t89 + t92;
t79 = sin(qJ(4));
t82 = cos(qJ(4));
t77 = sin(pkin(6));
t93 = t77 * t84;
t62 = t70 * t82 + t79 * t93;
t94 = t77 * t83;
t68 = t78 * t79 + t82 * t94;
t59 = atan2(t62, t68);
t56 = sin(t59);
t57 = cos(t59);
t50 = t56 * t62 + t57 * t68;
t49 = 0.1e1 / t50 ^ 2;
t90 = t84 * t80;
t91 = t81 * t83;
t85 = t78 * t91 + t90;
t95 = t77 * t81;
t60 = t79 * t95 - t85 * t82;
t102 = t49 * t60;
t61 = t85 * t79 + t82 * t95;
t72 = -t78 * t92 + t89;
t76 = pkin(11) + qJ(6);
t74 = sin(t76);
t75 = cos(t76);
t55 = t61 * t75 + t72 * t74;
t53 = 0.1e1 / t55 ^ 2;
t54 = t61 * t74 - t72 * t75;
t101 = t53 * t54;
t100 = t57 * t62;
t99 = t60 ^ 2 * t49;
t67 = 0.1e1 / t68 ^ 2;
t98 = t62 * t67;
t97 = t72 * t79;
t96 = t77 * t80;
t88 = t54 ^ 2 * t53 + 0.1e1;
t87 = -t56 * t68 + t100;
t86 = -t70 * t79 + t82 * t93;
t71 = t78 * t90 + t91;
t69 = t78 * t82 - t79 * t94;
t66 = 0.1e1 / t68;
t58 = 0.1e1 / (t62 ^ 2 * t67 + 0.1e1);
t52 = 0.1e1 / t55;
t51 = 0.1e1 / t88;
t48 = 0.1e1 / t50;
t47 = 0.1e1 / (0.1e1 + t99);
t46 = (t66 * t71 + t96 * t98) * t82 * t58;
t45 = (t66 * t86 - t69 * t98) * t58;
t1 = [-t60 * t66 * t58, t46, 0, t45, 0, 0; (t62 * t48 - (-t56 + (-t66 * t100 + t56) * t58) * t99) * t47 (-t72 * t82 * t48 - ((t56 * t71 - t57 * t96) * t82 + t87 * t46) * t102) * t47, 0 (t61 * t48 - (t87 * t45 + t56 * t86 + t57 * t69) * t102) * t47, 0, 0; ((t71 * t75 + t74 * t86) * t52 - (-t71 * t74 + t75 * t86) * t101) * t51 ((t74 * t97 + t85 * t75) * t52 - (-t85 * t74 + t75 * t97) * t101) * t51, 0 (t75 * t101 - t74 * t52) * t60 * t51, 0, t88 * t51;];
Ja_rot  = t1;
