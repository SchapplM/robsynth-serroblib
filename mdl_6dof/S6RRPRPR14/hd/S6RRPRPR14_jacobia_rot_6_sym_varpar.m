% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR14_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:45:26
% EndTime: 2019-02-26 21:45:26
% DurationCPUTime: 0.18s
% Computational Cost: add. (459->36), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->54)
t74 = sin(qJ(4));
t78 = cos(qJ(4));
t72 = cos(pkin(6));
t79 = cos(qJ(2));
t80 = cos(qJ(1));
t86 = t80 * t79;
t75 = sin(qJ(2));
t76 = sin(qJ(1));
t89 = t76 * t75;
t82 = -t72 * t86 + t89;
t71 = sin(pkin(6));
t90 = t71 * t80;
t61 = -t82 * t74 + t78 * t90;
t100 = t74 * t90 + t82 * t78;
t91 = t71 * t79;
t66 = t72 * t78 - t74 * t91;
t56 = atan2(t61, t66);
t53 = sin(t56);
t54 = cos(t56);
t47 = t53 * t61 + t54 * t66;
t46 = 0.1e1 / t47 ^ 2;
t87 = t80 * t75;
t88 = t76 * t79;
t68 = t72 * t88 + t87;
t92 = t71 * t76;
t58 = t68 * t74 + t78 * t92;
t99 = t46 * t58;
t57 = -t68 * t78 + t74 * t92;
t69 = -t72 * t89 + t86;
t73 = sin(qJ(6));
t77 = cos(qJ(6));
t52 = t57 * t73 + t69 * t77;
t50 = 0.1e1 / t52 ^ 2;
t51 = -t57 * t77 + t69 * t73;
t98 = t50 * t51;
t97 = t54 * t61;
t96 = t58 ^ 2 * t46;
t64 = 0.1e1 / t66 ^ 2;
t95 = t61 * t64;
t94 = t69 * t78;
t93 = t71 * t75;
t83 = t51 ^ 2 * t50 + 0.1e1;
t81 = -t53 * t66 + t97;
t67 = t72 * t87 + t88;
t65 = -t72 * t74 - t78 * t91;
t63 = 0.1e1 / t66;
t55 = 0.1e1 / (t61 ^ 2 * t64 + 0.1e1);
t49 = 0.1e1 / t52;
t48 = 0.1e1 / t83;
t45 = 0.1e1 / t47;
t44 = 0.1e1 / (0.1e1 + t96);
t43 = (-t63 * t67 - t93 * t95) * t74 * t55;
t42 = (-t100 * t63 - t65 * t95) * t55;
t1 = [-t58 * t63 * t55, t43, 0, t42, 0, 0; (t61 * t45 - (-t53 + (-t63 * t97 + t53) * t55) * t96) * t44 (t69 * t74 * t45 - ((-t53 * t67 + t54 * t93) * t74 + t81 * t43) * t99) * t44, 0 (-t57 * t45 - (-t100 * t53 + t81 * t42 + t54 * t65) * t99) * t44, 0, 0; ((-t100 * t77 - t67 * t73) * t49 - (t100 * t73 - t67 * t77) * t98) * t48 ((-t68 * t73 + t77 * t94) * t49 - (-t68 * t77 - t73 * t94) * t98) * t48, 0 (-t77 * t49 - t73 * t98) * t58 * t48, 0, t83 * t48;];
Ja_rot  = t1;
