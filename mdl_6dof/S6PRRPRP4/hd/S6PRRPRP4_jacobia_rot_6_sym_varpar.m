% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRP4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:03:04
% EndTime: 2019-02-26 20:03:04
% DurationCPUTime: 0.11s
% Computational Cost: add. (334->30), mult. (949->78), div. (65->9), fcn. (1349->13), ass. (0->50)
t72 = sin(pkin(10));
t74 = cos(pkin(10));
t81 = cos(qJ(2));
t75 = cos(pkin(6));
t78 = sin(qJ(2));
t85 = t75 * t78;
t66 = t72 * t81 + t74 * t85;
t80 = cos(qJ(3));
t73 = sin(pkin(6));
t77 = sin(qJ(3));
t88 = t73 * t77;
t59 = t66 * t80 - t74 * t88;
t87 = t73 * t80;
t70 = t75 * t77 + t78 * t87;
t57 = atan2(-t59, t70);
t54 = sin(t57);
t55 = cos(t57);
t48 = -t54 * t59 + t55 * t70;
t47 = 0.1e1 / t48 ^ 2;
t68 = -t72 * t85 + t74 * t81;
t62 = t68 * t80 + t72 * t88;
t93 = t47 * t62;
t61 = t68 * t77 - t72 * t87;
t76 = sin(qJ(5));
t84 = t75 * t81;
t67 = t72 * t84 + t74 * t78;
t79 = cos(qJ(5));
t89 = t67 * t79;
t53 = t61 * t76 + t89;
t51 = 0.1e1 / t53 ^ 2;
t90 = t67 * t76;
t52 = -t61 * t79 + t90;
t92 = t51 * t52;
t64 = 0.1e1 / t70 ^ 2;
t91 = t59 * t64;
t86 = t73 * t81;
t83 = t52 ^ 2 * t51 + 0.1e1;
t82 = -t54 * t70 - t55 * t59;
t69 = t75 * t80 - t78 * t88;
t65 = -t72 * t78 + t74 * t84;
t63 = 0.1e1 / t70;
t58 = t66 * t77 + t74 * t87;
t56 = 0.1e1 / (t59 ^ 2 * t64 + 0.1e1);
t50 = 0.1e1 / t53;
t49 = 0.1e1 / t83;
t46 = 0.1e1 / t48;
t45 = 0.1e1 / (t62 ^ 2 * t47 + 0.1e1);
t44 = (-t63 * t65 + t86 * t91) * t80 * t56;
t43 = (t58 * t63 + t69 * t91) * t56;
t1 = [0, t44, t43, 0, 0, 0; 0 (-t67 * t80 * t46 - ((-t54 * t65 + t55 * t86) * t80 + t82 * t44) * t93) * t45 (-t61 * t46 - (t82 * t43 + t54 * t58 + t55 * t69) * t93) * t45, 0, 0, 0; 0 ((t68 * t76 + t77 * t89) * t50 - (t68 * t79 - t77 * t90) * t92) * t49 (-t50 * t79 - t76 * t92) * t62 * t49, 0, t83 * t49, 0;];
Ja_rot  = t1;
