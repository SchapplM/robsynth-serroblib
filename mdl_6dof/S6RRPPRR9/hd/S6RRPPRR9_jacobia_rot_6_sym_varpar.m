% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:12
% EndTime: 2019-02-26 21:33:13
% DurationCPUTime: 0.17s
% Computational Cost: add. (459->36), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->55)
t75 = cos(pkin(6));
t79 = sin(qJ(1));
t82 = cos(qJ(2));
t89 = t79 * t82;
t78 = sin(qJ(2));
t83 = cos(qJ(1));
t90 = t78 * t83;
t70 = t75 * t90 + t89;
t77 = sin(qJ(5));
t81 = cos(qJ(5));
t74 = sin(pkin(6));
t92 = t74 * t83;
t61 = t70 * t81 + t77 * t92;
t94 = t74 * t81;
t67 = t75 * t77 - t78 * t94;
t58 = atan2(t61, t67);
t55 = sin(t58);
t56 = cos(t58);
t49 = t55 * t61 + t56 * t67;
t48 = 0.1e1 / t49 ^ 2;
t88 = t82 * t83;
t91 = t78 * t79;
t86 = -t75 * t91 + t88;
t95 = t74 * t77;
t59 = t79 * t95 - t86 * t81;
t102 = t48 * t59;
t101 = t48 * t59 ^ 2;
t60 = t86 * t77 + t79 * t94;
t80 = cos(qJ(6));
t71 = -t75 * t89 - t90;
t76 = sin(qJ(6));
t97 = t71 * t76;
t54 = t60 * t80 + t97;
t52 = 0.1e1 / t54 ^ 2;
t96 = t71 * t80;
t53 = t60 * t76 - t96;
t100 = t52 * t53;
t99 = t56 * t61;
t66 = 0.1e1 / t67 ^ 2;
t98 = t61 * t66;
t93 = t74 * t82;
t87 = t52 * t53 ^ 2 + 0.1e1;
t85 = -t55 * t67 + t99;
t84 = -t70 * t77 + t81 * t92;
t69 = -t75 * t88 + t91;
t68 = t75 * t81 + t78 * t95;
t65 = 0.1e1 / t67;
t57 = 0.1e1 / (t61 ^ 2 * t66 + 0.1e1);
t51 = 0.1e1 / t54;
t50 = 0.1e1 / t87;
t47 = 0.1e1 / t49;
t46 = 0.1e1 / (0.1e1 + t101);
t45 = (-t65 * t69 + t93 * t98) * t81 * t57;
t44 = (t65 * t84 - t68 * t98) * t57;
t1 = [-t59 * t65 * t57, t45, 0, 0, t44, 0; (t61 * t47 - (-t55 + (-t65 * t99 + t55) * t57) * t101) * t46 (-t71 * t81 * t47 - ((-t55 * t69 - t56 * t93) * t81 + t85 * t45) * t102) * t46, 0, 0 (t60 * t47 - (t85 * t44 + t55 * t84 + t56 * t68) * t102) * t46, 0; ((-t69 * t80 + t76 * t84) * t51 - (t69 * t76 + t80 * t84) * t100) * t50 ((t77 * t97 + t86 * t80) * t51 - (-t86 * t76 + t77 * t96) * t100) * t50, 0, 0 (t80 * t100 - t76 * t51) * t59 * t50, t87 * t50;];
Ja_rot  = t1;
