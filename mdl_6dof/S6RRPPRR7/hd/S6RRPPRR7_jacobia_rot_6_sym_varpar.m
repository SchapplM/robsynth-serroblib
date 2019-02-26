% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR7
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
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:57
% EndTime: 2019-02-26 21:31:57
% DurationCPUTime: 0.18s
% Computational Cost: add. (459->37), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->54)
t73 = cos(pkin(6));
t80 = cos(qJ(2));
t81 = cos(qJ(1));
t84 = t81 * t80;
t76 = sin(qJ(2));
t77 = sin(qJ(1));
t87 = t77 * t76;
t66 = -t73 * t84 + t87;
t75 = sin(qJ(5));
t79 = cos(qJ(5));
t72 = sin(pkin(6));
t88 = t72 * t81;
t56 = t66 * t75 - t79 * t88;
t89 = t72 * t80;
t64 = t73 * t79 - t75 * t89;
t55 = atan2(-t56, t64);
t52 = sin(t55);
t53 = cos(t55);
t46 = -t52 * t56 + t53 * t64;
t45 = 0.1e1 / t46 ^ 2;
t85 = t81 * t76;
t86 = t77 * t80;
t68 = t73 * t86 + t85;
t90 = t72 * t77;
t60 = t68 * t75 + t79 * t90;
t97 = t45 * t60;
t61 = t68 * t79 - t75 * t90;
t69 = -t73 * t87 + t84;
t74 = sin(qJ(6));
t78 = cos(qJ(6));
t51 = t61 * t78 + t69 * t74;
t49 = 0.1e1 / t51 ^ 2;
t50 = t61 * t74 - t69 * t78;
t96 = t49 * t50;
t95 = t53 * t56;
t63 = 0.1e1 / t64 ^ 2;
t94 = t56 * t63;
t93 = t60 ^ 2 * t45;
t92 = t69 * t79;
t91 = t72 * t76;
t83 = t50 ^ 2 * t49 + 0.1e1;
t82 = -t52 * t64 - t95;
t58 = t66 * t79 + t75 * t88;
t67 = t73 * t85 + t86;
t65 = -t73 * t75 - t79 * t89;
t62 = 0.1e1 / t64;
t54 = 0.1e1 / (t56 ^ 2 * t63 + 0.1e1);
t48 = 0.1e1 / t51;
t47 = 0.1e1 / t83;
t44 = 0.1e1 / t46;
t43 = 0.1e1 / (0.1e1 + t93);
t42 = (-t62 * t67 + t91 * t94) * t75 * t54;
t41 = (-t58 * t62 + t65 * t94) * t54;
t1 = [-t60 * t62 * t54, t42, 0, 0, t41, 0; (-t56 * t44 - (-t52 + (t62 * t95 + t52) * t54) * t93) * t43 (t69 * t75 * t44 - ((-t52 * t67 + t53 * t91) * t75 + t82 * t42) * t97) * t43, 0, 0 (t61 * t44 - (t82 * t41 - t52 * t58 + t53 * t65) * t97) * t43, 0; ((-t58 * t74 + t67 * t78) * t48 - (-t58 * t78 - t67 * t74) * t96) * t47 ((t68 * t78 + t74 * t92) * t48 - (-t68 * t74 + t78 * t92) * t96) * t47, 0, 0 (-t74 * t48 + t78 * t96) * t60 * t47, t83 * t47;];
Ja_rot  = t1;
