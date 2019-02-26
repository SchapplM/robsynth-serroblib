% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:06:23
% EndTime: 2019-02-26 22:06:23
% DurationCPUTime: 0.22s
% Computational Cost: add. (878->38), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
t81 = cos(pkin(6));
t84 = sin(qJ(1));
t86 = cos(qJ(2));
t91 = t84 * t86;
t83 = sin(qJ(2));
t87 = cos(qJ(1));
t93 = t83 * t87;
t73 = t81 * t93 + t91;
t79 = qJ(3) + pkin(11);
t77 = sin(t79);
t78 = cos(t79);
t80 = sin(pkin(6));
t94 = t80 * t87;
t63 = t73 * t78 - t77 * t94;
t97 = t80 * t83;
t71 = t77 * t81 + t78 * t97;
t61 = atan2(-t63, t71);
t54 = sin(t61);
t55 = cos(t61);
t52 = -t54 * t63 + t55 * t71;
t51 = 0.1e1 / t52 ^ 2;
t90 = t86 * t87;
t92 = t84 * t83;
t75 = -t81 * t92 + t90;
t96 = t80 * t84;
t67 = t75 * t78 + t77 * t96;
t104 = t51 * t67;
t103 = t51 * t67 ^ 2;
t102 = t55 * t63;
t66 = t75 * t77 - t78 * t96;
t82 = sin(qJ(6));
t74 = t81 * t91 + t93;
t85 = cos(qJ(6));
t98 = t74 * t85;
t60 = t66 * t82 + t98;
t57 = 0.1e1 / t60 ^ 2;
t99 = t74 * t82;
t59 = -t66 * t85 + t99;
t101 = t57 * t59;
t69 = 0.1e1 / t71 ^ 2;
t100 = t63 * t69;
t95 = t80 * t86;
t89 = t57 * t59 ^ 2 + 0.1e1;
t88 = -t54 * t71 - t102;
t62 = t73 * t77 + t78 * t94;
t72 = t81 * t90 - t92;
t70 = -t77 * t97 + t78 * t81;
t68 = 0.1e1 / t71;
t58 = 0.1e1 / (t63 ^ 2 * t69 + 0.1e1);
t56 = 0.1e1 / t60;
t53 = 0.1e1 / t89;
t50 = 0.1e1 / t52;
t49 = 0.1e1 / (0.1e1 + t103);
t48 = (t95 * t100 - t68 * t72) * t78 * t58;
t47 = (t70 * t100 + t62 * t68) * t58;
t1 = [-t67 * t68 * t58, t48, t47, 0, 0, 0; (-t63 * t50 - (-t54 + (t68 * t102 + t54) * t58) * t103) * t49 (-t74 * t78 * t50 - ((-t54 * t72 + t55 * t95) * t78 + t88 * t48) * t104) * t49 (-t66 * t50 - (t88 * t47 + t54 * t62 + t55 * t70) * t104) * t49, 0, 0, 0; ((t62 * t85 + t72 * t82) * t56 - (-t62 * t82 + t72 * t85) * t101) * t53 ((t75 * t82 + t77 * t98) * t56 - (t75 * t85 - t77 * t99) * t101) * t53 (-t82 * t101 - t56 * t85) * t67 * t53, 0, 0, t89 * t53;];
Ja_rot  = t1;
