% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRP3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:21
% EndTime: 2019-02-26 20:16:21
% DurationCPUTime: 0.12s
% Computational Cost: add. (427->31), mult. (1032->79), div. (70->9), fcn. (1461->13), ass. (0->51)
t80 = sin(pkin(11));
t82 = cos(pkin(11));
t87 = cos(qJ(2));
t83 = cos(pkin(6));
t85 = sin(qJ(2));
t91 = t83 * t85;
t71 = t80 * t87 + t82 * t91;
t84 = sin(qJ(3));
t81 = sin(pkin(6));
t86 = cos(qJ(3));
t93 = t81 * t86;
t63 = t71 * t84 + t82 * t93;
t94 = t81 * t84;
t74 = -t83 * t86 + t85 * t94;
t62 = atan2(-t63, t74);
t59 = sin(t62);
t60 = cos(t62);
t53 = -t59 * t63 + t60 * t74;
t52 = 0.1e1 / t53 ^ 2;
t73 = -t80 * t91 + t82 * t87;
t66 = t73 * t84 - t80 * t93;
t98 = t52 * t66;
t67 = t73 * t86 + t80 * t94;
t90 = t83 * t87;
t72 = t80 * t90 + t82 * t85;
t79 = qJ(4) + qJ(5);
t77 = sin(t79);
t78 = cos(t79);
t58 = t67 * t78 + t72 * t77;
t56 = 0.1e1 / t58 ^ 2;
t57 = t67 * t77 - t72 * t78;
t97 = t56 * t57;
t69 = 0.1e1 / t74 ^ 2;
t96 = t63 * t69;
t95 = t72 * t86;
t92 = t81 * t87;
t89 = t57 ^ 2 * t56 + 0.1e1;
t88 = -t59 * t74 - t60 * t63;
t75 = t83 * t84 + t85 * t93;
t70 = -t80 * t85 + t82 * t90;
t68 = 0.1e1 / t74;
t65 = t71 * t86 - t82 * t94;
t61 = 0.1e1 / (t63 ^ 2 * t69 + 0.1e1);
t55 = 0.1e1 / t58;
t54 = 0.1e1 / t89;
t51 = 0.1e1 / t53;
t50 = 0.1e1 / (t66 ^ 2 * t52 + 0.1e1);
t49 = (-t68 * t70 + t92 * t96) * t84 * t61;
t48 = (-t65 * t68 + t75 * t96) * t61;
t47 = t89 * t54;
t1 = [0, t49, t48, 0, 0, 0; 0 (-t72 * t84 * t51 - ((-t59 * t70 + t60 * t92) * t84 + t88 * t49) * t98) * t50 (t67 * t51 - (t88 * t48 - t59 * t65 + t60 * t75) * t98) * t50, 0, 0, 0; 0 ((-t73 * t78 - t77 * t95) * t55 - (t73 * t77 - t78 * t95) * t97) * t54 (-t55 * t77 + t78 * t97) * t66 * t54, t47, t47, 0;];
Ja_rot  = t1;
