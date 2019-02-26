% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR14_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:44
% EndTime: 2019-02-26 22:23:45
% DurationCPUTime: 0.17s
% Computational Cost: add. (459->37), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->55)
t72 = cos(pkin(6));
t75 = sin(qJ(2));
t80 = cos(qJ(1));
t84 = t80 * t75;
t76 = sin(qJ(1));
t79 = cos(qJ(2));
t85 = t76 * t79;
t67 = t72 * t84 + t85;
t74 = sin(qJ(3));
t78 = cos(qJ(3));
t71 = sin(pkin(6));
t87 = t71 * t80;
t57 = t67 * t78 - t74 * t87;
t89 = t71 * t78;
t65 = t72 * t74 + t75 * t89;
t55 = atan2(-t57, t65);
t52 = sin(t55);
t53 = cos(t55);
t46 = -t52 * t57 + t53 * t65;
t45 = 0.1e1 / t46 ^ 2;
t83 = t80 * t79;
t86 = t76 * t75;
t69 = -t72 * t86 + t83;
t90 = t71 * t74;
t61 = t69 * t78 + t76 * t90;
t97 = t45 * t61;
t60 = t69 * t74 - t76 * t89;
t73 = sin(qJ(5));
t68 = t72 * t85 + t84;
t77 = cos(qJ(5));
t91 = t68 * t77;
t51 = t60 * t73 + t91;
t49 = 0.1e1 / t51 ^ 2;
t92 = t68 * t73;
t50 = -t60 * t77 + t92;
t96 = t49 * t50;
t95 = t53 * t57;
t63 = 0.1e1 / t65 ^ 2;
t94 = t57 * t63;
t93 = t61 ^ 2 * t45;
t88 = t71 * t79;
t82 = t50 ^ 2 * t49 + 0.1e1;
t81 = -t52 * t65 - t95;
t56 = t67 * t74 + t78 * t87;
t66 = t72 * t83 - t86;
t64 = t72 * t78 - t75 * t90;
t62 = 0.1e1 / t65;
t54 = 0.1e1 / (t57 ^ 2 * t63 + 0.1e1);
t48 = 0.1e1 / t51;
t47 = 0.1e1 / t82;
t44 = 0.1e1 / t46;
t43 = 0.1e1 / (0.1e1 + t93);
t42 = (-t62 * t66 + t88 * t94) * t78 * t54;
t41 = (t56 * t62 + t64 * t94) * t54;
t1 = [-t61 * t62 * t54, t42, t41, 0, 0, 0; (-t57 * t44 - (-t52 + (t62 * t95 + t52) * t54) * t93) * t43 (-t68 * t78 * t44 - ((-t52 * t66 + t53 * t88) * t78 + t81 * t42) * t97) * t43 (-t60 * t44 - (t81 * t41 + t52 * t56 + t53 * t64) * t97) * t43, 0, 0, 0; ((t56 * t77 + t66 * t73) * t48 - (-t56 * t73 + t66 * t77) * t96) * t47 ((t69 * t73 + t74 * t91) * t48 - (t69 * t77 - t74 * t92) * t96) * t47 (-t77 * t48 - t73 * t96) * t61 * t47, 0, t82 * t47, 0;];
Ja_rot  = t1;
