% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP11_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:14:59
% EndTime: 2019-02-26 22:14:59
% DurationCPUTime: 0.16s
% Computational Cost: add. (459->37), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->55)
t78 = cos(pkin(6));
t81 = sin(qJ(2));
t86 = cos(qJ(1));
t90 = t86 * t81;
t82 = sin(qJ(1));
t85 = cos(qJ(2));
t91 = t82 * t85;
t73 = t78 * t90 + t91;
t80 = sin(qJ(3));
t84 = cos(qJ(3));
t77 = sin(pkin(6));
t93 = t77 * t86;
t63 = t73 * t84 - t80 * t93;
t95 = t77 * t84;
t71 = t78 * t80 + t81 * t95;
t61 = atan2(-t63, t71);
t58 = sin(t61);
t59 = cos(t61);
t52 = -t58 * t63 + t59 * t71;
t51 = 0.1e1 / t52 ^ 2;
t89 = t86 * t85;
t92 = t82 * t81;
t75 = -t78 * t92 + t89;
t96 = t77 * t80;
t67 = t75 * t84 + t82 * t96;
t103 = t51 * t67;
t66 = t75 * t80 - t82 * t95;
t79 = sin(qJ(5));
t74 = t78 * t91 + t90;
t83 = cos(qJ(5));
t97 = t74 * t83;
t57 = t66 * t79 + t97;
t55 = 0.1e1 / t57 ^ 2;
t98 = t74 * t79;
t56 = -t66 * t83 + t98;
t102 = t55 * t56;
t101 = t59 * t63;
t69 = 0.1e1 / t71 ^ 2;
t100 = t63 * t69;
t99 = t67 ^ 2 * t51;
t94 = t77 * t85;
t88 = t56 ^ 2 * t55 + 0.1e1;
t87 = -t58 * t71 - t101;
t62 = t73 * t80 + t84 * t93;
t72 = t78 * t89 - t92;
t70 = t78 * t84 - t81 * t96;
t68 = 0.1e1 / t71;
t60 = 0.1e1 / (t63 ^ 2 * t69 + 0.1e1);
t54 = 0.1e1 / t57;
t53 = 0.1e1 / t88;
t50 = 0.1e1 / t52;
t49 = 0.1e1 / (0.1e1 + t99);
t48 = (t94 * t100 - t68 * t72) * t84 * t60;
t47 = (t70 * t100 + t62 * t68) * t60;
t1 = [-t67 * t68 * t60, t48, t47, 0, 0, 0; (-t63 * t50 - (-t58 + (t68 * t101 + t58) * t60) * t99) * t49 (-t74 * t84 * t50 - ((-t58 * t72 + t59 * t94) * t84 + t87 * t48) * t103) * t49 (-t66 * t50 - (t87 * t47 + t58 * t62 + t59 * t70) * t103) * t49, 0, 0, 0; ((t62 * t83 + t72 * t79) * t54 - (-t62 * t79 + t72 * t83) * t102) * t53 ((t75 * t79 + t80 * t97) * t54 - (t75 * t83 - t80 * t98) * t102) * t53 (-t79 * t102 - t83 * t54) * t67 * t53, 0, t88 * t53, 0;];
Ja_rot  = t1;
