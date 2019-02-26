% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:39
% EndTime: 2019-02-26 22:07:40
% DurationCPUTime: 0.18s
% Computational Cost: add. (459->37), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->55)
t76 = cos(pkin(6));
t79 = sin(qJ(2));
t84 = cos(qJ(1));
t88 = t84 * t79;
t80 = sin(qJ(1));
t83 = cos(qJ(2));
t89 = t80 * t83;
t69 = t76 * t88 + t89;
t78 = sin(qJ(3));
t82 = cos(qJ(3));
t75 = sin(pkin(6));
t91 = t75 * t84;
t59 = t69 * t82 - t78 * t91;
t93 = t75 * t82;
t67 = t76 * t78 + t79 * t93;
t57 = atan2(-t59, t67);
t54 = sin(t57);
t55 = cos(t57);
t48 = -t54 * t59 + t55 * t67;
t47 = 0.1e1 / t48 ^ 2;
t87 = t84 * t83;
t90 = t80 * t79;
t71 = -t76 * t90 + t87;
t94 = t75 * t78;
t63 = t71 * t82 + t80 * t94;
t101 = t47 * t63;
t62 = t71 * t78 - t80 * t93;
t81 = cos(qJ(6));
t70 = -t76 * t89 - t88;
t77 = sin(qJ(6));
t96 = t70 * t77;
t53 = t62 * t81 + t96;
t51 = 0.1e1 / t53 ^ 2;
t95 = t70 * t81;
t52 = t62 * t77 - t95;
t100 = t51 * t52;
t99 = t55 * t59;
t65 = 0.1e1 / t67 ^ 2;
t98 = t59 * t65;
t97 = t63 ^ 2 * t47;
t92 = t75 * t83;
t86 = t52 ^ 2 * t51 + 0.1e1;
t85 = -t54 * t67 - t99;
t58 = t69 * t78 + t82 * t91;
t68 = -t76 * t87 + t90;
t66 = t76 * t82 - t79 * t94;
t64 = 0.1e1 / t67;
t56 = 0.1e1 / (t59 ^ 2 * t65 + 0.1e1);
t50 = 0.1e1 / t53;
t49 = 0.1e1 / t86;
t46 = 0.1e1 / t48;
t45 = 0.1e1 / (0.1e1 + t97);
t44 = (t64 * t68 + t92 * t98) * t82 * t56;
t43 = (t58 * t64 + t66 * t98) * t56;
t1 = [-t63 * t64 * t56, t44, t43, 0, 0, 0; (-t59 * t46 - (-t54 + (t64 * t99 + t54) * t56) * t97) * t45 (t70 * t82 * t46 - ((t54 * t68 + t55 * t92) * t82 + t85 * t44) * t101) * t45 (-t62 * t46 - (t43 * t85 + t54 * t58 + t55 * t66) * t101) * t45, 0, 0, 0; ((-t58 * t77 - t68 * t81) * t50 - (-t58 * t81 + t68 * t77) * t100) * t49 ((t71 * t81 + t78 * t96) * t50 - (-t71 * t77 + t78 * t95) * t100) * t49 (-t81 * t100 + t77 * t50) * t63 * t49, 0, 0, t86 * t49;];
Ja_rot  = t1;
