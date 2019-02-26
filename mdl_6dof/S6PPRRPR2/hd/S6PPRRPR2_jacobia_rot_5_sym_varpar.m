% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRPR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRPR2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:56
% EndTime: 2019-02-26 19:40:57
% DurationCPUTime: 0.24s
% Computational Cost: add. (826->40), mult. (2382->97), div. (57->9), fcn. (3277->15), ass. (0->58)
t81 = sin(pkin(12));
t82 = sin(pkin(11));
t85 = cos(pkin(12));
t86 = cos(pkin(11));
t88 = cos(pkin(6));
t99 = t86 * t88;
t77 = t81 * t99 + t82 * t85;
t90 = sin(qJ(3));
t92 = cos(qJ(3));
t83 = sin(pkin(7));
t84 = sin(pkin(6));
t101 = t84 * t83;
t87 = cos(pkin(7));
t95 = -t81 * t82 + t85 * t99;
t93 = -t101 * t86 + t87 * t95;
t64 = t77 * t92 + t90 * t93;
t89 = sin(qJ(4));
t91 = cos(qJ(4));
t100 = t84 * t87;
t94 = -t100 * t86 - t83 * t95;
t56 = t64 * t89 - t91 * t94;
t102 = t83 * t88;
t73 = t90 * t102 + (t85 * t87 * t90 + t81 * t92) * t84;
t76 = -t101 * t85 + t87 * t88;
t69 = t73 * t89 - t76 * t91;
t55 = atan2(-t56, t69);
t52 = sin(t55);
t53 = cos(t55);
t50 = -t52 * t56 + t53 * t69;
t49 = 0.1e1 / t50 ^ 2;
t103 = t82 * t88;
t78 = -t103 * t85 - t81 * t86;
t79 = -t103 * t81 + t85 * t86;
t97 = t82 * t101;
t66 = t79 * t92 + (t78 * t87 + t97) * t90;
t74 = t100 * t82 - t78 * t83;
t59 = t66 * t89 - t74 * t91;
t105 = t49 * t59;
t68 = 0.1e1 / t69 ^ 2;
t104 = t56 * t68;
t98 = t87 * t92;
t96 = -t52 * t69 - t53 * t56;
t72 = t92 * t102 + (-t81 * t90 + t85 * t98) * t84;
t70 = t73 * t91 + t76 * t89;
t67 = 0.1e1 / t69;
t65 = -t78 * t98 + t79 * t90 - t92 * t97;
t63 = -t77 * t90 + t92 * t93;
t62 = 0.1e1 / t65 ^ 2;
t61 = 0.1e1 / t65;
t60 = t66 * t91 + t74 * t89;
t58 = t64 * t91 + t89 * t94;
t54 = 0.1e1 / (t56 ^ 2 * t68 + 0.1e1);
t51 = 0.1e1 / (t60 ^ 2 * t62 + 0.1e1);
t48 = 0.1e1 / t50;
t47 = 0.1e1 / (t49 * t59 ^ 2 + 0.1e1);
t46 = (t104 * t72 - t63 * t67) * t89 * t54;
t45 = (t104 * t70 - t58 * t67) * t54;
t1 = [0, 0, t46, t45, 0, 0; 0, 0 (-t65 * t89 * t48 - ((-t52 * t63 + t53 * t72) * t89 + t96 * t46) * t105) * t47 (t60 * t48 - (t45 * t96 - t52 * t58 + t53 * t70) * t105) * t47, 0, 0; 0, 0 (-t60 * t62 * t66 - t61 * t65 * t91) * t51, -t59 * t61 * t51, 0, 0;];
Ja_rot  = t1;
