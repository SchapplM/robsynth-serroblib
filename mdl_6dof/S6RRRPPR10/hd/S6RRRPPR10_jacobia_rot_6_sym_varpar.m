% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR10_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:52
% EndTime: 2019-02-26 22:08:52
% DurationCPUTime: 0.16s
% Computational Cost: add. (525->38), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->55)
t79 = cos(pkin(6));
t81 = sin(qJ(2));
t85 = cos(qJ(1));
t89 = t85 * t81;
t82 = sin(qJ(1));
t84 = cos(qJ(2));
t90 = t82 * t84;
t71 = t79 * t89 + t90;
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t78 = sin(pkin(6));
t92 = t78 * t85;
t61 = t71 * t83 - t80 * t92;
t94 = t78 * t83;
t69 = t79 * t80 + t81 * t94;
t59 = atan2(-t61, t69);
t56 = sin(t59);
t57 = cos(t59);
t50 = -t56 * t61 + t57 * t69;
t49 = 0.1e1 / t50 ^ 2;
t88 = t85 * t84;
t91 = t82 * t81;
t73 = -t79 * t91 + t88;
t95 = t78 * t80;
t65 = t73 * t83 + t82 * t95;
t101 = t49 * t65;
t64 = t73 * t80 - t82 * t94;
t72 = t79 * t90 + t89;
t77 = pkin(11) + qJ(6);
t75 = sin(t77);
t76 = cos(t77);
t55 = t64 * t75 + t72 * t76;
t53 = 0.1e1 / t55 ^ 2;
t54 = -t64 * t76 + t72 * t75;
t100 = t53 * t54;
t99 = t57 * t61;
t67 = 0.1e1 / t69 ^ 2;
t98 = t61 * t67;
t97 = t65 ^ 2 * t49;
t96 = t72 * t80;
t93 = t78 * t84;
t87 = t54 ^ 2 * t53 + 0.1e1;
t86 = -t56 * t69 - t99;
t60 = t71 * t80 + t83 * t92;
t70 = t79 * t88 - t91;
t68 = t79 * t83 - t81 * t95;
t66 = 0.1e1 / t69;
t58 = 0.1e1 / (t61 ^ 2 * t67 + 0.1e1);
t52 = 0.1e1 / t55;
t51 = 0.1e1 / t87;
t48 = 0.1e1 / t50;
t47 = 0.1e1 / (0.1e1 + t97);
t46 = (-t66 * t70 + t93 * t98) * t83 * t58;
t45 = (t60 * t66 + t68 * t98) * t58;
t1 = [-t65 * t66 * t58, t46, t45, 0, 0, 0; (-t61 * t48 - (-t56 + (t66 * t99 + t56) * t58) * t97) * t47 (-t72 * t83 * t48 - ((-t56 * t70 + t57 * t93) * t83 + t86 * t46) * t101) * t47 (-t64 * t48 - (t86 * t45 + t56 * t60 + t57 * t68) * t101) * t47, 0, 0, 0; ((t60 * t76 + t70 * t75) * t52 - (-t60 * t75 + t70 * t76) * t100) * t51 ((t73 * t75 + t76 * t96) * t52 - (t73 * t76 - t75 * t96) * t100) * t51 (-t75 * t100 - t76 * t52) * t65 * t51, 0, 0, t87 * t51;];
Ja_rot  = t1;
