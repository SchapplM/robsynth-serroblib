% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:06:27
% EndTime: 2019-02-26 20:06:27
% DurationCPUTime: 0.12s
% Computational Cost: add. (489->31), mult. (1032->79), div. (70->9), fcn. (1461->13), ass. (0->51)
t84 = sin(pkin(11));
t86 = cos(pkin(11));
t91 = cos(qJ(2));
t87 = cos(pkin(6));
t89 = sin(qJ(2));
t95 = t87 * t89;
t75 = t84 * t91 + t86 * t95;
t88 = sin(qJ(3));
t85 = sin(pkin(6));
t90 = cos(qJ(3));
t97 = t85 * t90;
t67 = t75 * t88 + t86 * t97;
t98 = t85 * t88;
t78 = -t87 * t90 + t89 * t98;
t66 = atan2(-t67, t78);
t63 = sin(t66);
t64 = cos(t66);
t57 = -t63 * t67 + t64 * t78;
t56 = 0.1e1 / t57 ^ 2;
t77 = -t84 * t95 + t86 * t91;
t70 = t77 * t88 - t84 * t97;
t102 = t56 * t70;
t71 = t77 * t90 + t84 * t98;
t94 = t87 * t91;
t76 = t84 * t94 + t86 * t89;
t83 = pkin(12) + qJ(5) + qJ(6);
t81 = sin(t83);
t82 = cos(t83);
t62 = t71 * t82 + t76 * t81;
t60 = 0.1e1 / t62 ^ 2;
t61 = t71 * t81 - t76 * t82;
t101 = t60 * t61;
t73 = 0.1e1 / t78 ^ 2;
t100 = t67 * t73;
t99 = t76 * t90;
t96 = t85 * t91;
t93 = t61 ^ 2 * t60 + 0.1e1;
t92 = -t63 * t78 - t64 * t67;
t79 = t87 * t88 + t89 * t97;
t74 = -t84 * t89 + t86 * t94;
t72 = 0.1e1 / t78;
t69 = t75 * t90 - t86 * t98;
t65 = 0.1e1 / (t67 ^ 2 * t73 + 0.1e1);
t59 = 0.1e1 / t62;
t58 = 0.1e1 / t93;
t55 = 0.1e1 / t57;
t54 = 0.1e1 / (t70 ^ 2 * t56 + 0.1e1);
t53 = (t96 * t100 - t72 * t74) * t88 * t65;
t52 = (t79 * t100 - t69 * t72) * t65;
t51 = t93 * t58;
t1 = [0, t53, t52, 0, 0, 0; 0 (-t76 * t88 * t55 - ((-t63 * t74 + t64 * t96) * t88 + t92 * t53) * t102) * t54 (t71 * t55 - (t92 * t52 - t63 * t69 + t64 * t79) * t102) * t54, 0, 0, 0; 0 ((-t77 * t82 - t81 * t99) * t59 - (t77 * t81 - t82 * t99) * t101) * t58 (t82 * t101 - t59 * t81) * t70 * t58, 0, t51, t51;];
Ja_rot  = t1;
