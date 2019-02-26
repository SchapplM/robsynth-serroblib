% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:37
% EndTime: 2019-02-26 20:07:37
% DurationCPUTime: 0.12s
% Computational Cost: add. (427->31), mult. (1032->79), div. (70->9), fcn. (1461->13), ass. (0->51)
t82 = sin(pkin(11));
t84 = cos(pkin(11));
t89 = cos(qJ(2));
t85 = cos(pkin(6));
t87 = sin(qJ(2));
t93 = t85 * t87;
t73 = t82 * t89 + t84 * t93;
t88 = cos(qJ(3));
t83 = sin(pkin(6));
t86 = sin(qJ(3));
t96 = t83 * t86;
t66 = t73 * t88 - t84 * t96;
t95 = t83 * t88;
t77 = t85 * t86 + t87 * t95;
t64 = atan2(-t66, t77);
t61 = sin(t64);
t62 = cos(t64);
t55 = -t61 * t66 + t62 * t77;
t54 = 0.1e1 / t55 ^ 2;
t75 = -t82 * t93 + t84 * t89;
t69 = t75 * t88 + t82 * t96;
t100 = t54 * t69;
t68 = t75 * t86 - t82 * t95;
t92 = t85 * t89;
t74 = t82 * t92 + t84 * t87;
t81 = qJ(5) + qJ(6);
t79 = sin(t81);
t80 = cos(t81);
t60 = t68 * t79 + t74 * t80;
t58 = 0.1e1 / t60 ^ 2;
t59 = -t68 * t80 + t74 * t79;
t99 = t58 * t59;
t71 = 0.1e1 / t77 ^ 2;
t98 = t66 * t71;
t97 = t74 * t86;
t94 = t83 * t89;
t91 = t59 ^ 2 * t58 + 0.1e1;
t90 = -t61 * t77 - t62 * t66;
t76 = t85 * t88 - t87 * t96;
t72 = -t82 * t87 + t84 * t92;
t70 = 0.1e1 / t77;
t65 = t73 * t86 + t84 * t95;
t63 = 0.1e1 / (t66 ^ 2 * t71 + 0.1e1);
t57 = 0.1e1 / t60;
t56 = 0.1e1 / t91;
t53 = 0.1e1 / t55;
t52 = 0.1e1 / (t69 ^ 2 * t54 + 0.1e1);
t51 = (-t70 * t72 + t94 * t98) * t88 * t63;
t50 = (t65 * t70 + t76 * t98) * t63;
t49 = t91 * t56;
t1 = [0, t51, t50, 0, 0, 0; 0 (-t74 * t88 * t53 - ((-t61 * t72 + t62 * t94) * t88 + t90 * t51) * t100) * t52 (-t68 * t53 - (t90 * t50 + t61 * t65 + t62 * t76) * t100) * t52, 0, 0, 0; 0 ((t75 * t79 + t80 * t97) * t57 - (t75 * t80 - t79 * t97) * t99) * t56 (-t57 * t80 - t79 * t99) * t69 * t56, 0, t49, t49;];
Ja_rot  = t1;
