% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:55
% EndTime: 2019-02-26 19:56:55
% DurationCPUTime: 0.13s
% Computational Cost: add. (427->30), mult. (1032->80), div. (70->9), fcn. (1461->13), ass. (0->52)
t79 = sin(pkin(11));
t81 = cos(pkin(11));
t84 = sin(qJ(2));
t82 = cos(pkin(6));
t86 = cos(qJ(2));
t89 = t82 * t86;
t70 = t79 * t84 - t81 * t89;
t85 = cos(qJ(4));
t80 = sin(pkin(6));
t83 = sin(qJ(4));
t94 = t80 * t83;
t65 = t70 * t85 + t81 * t94;
t91 = t80 * t86;
t74 = t82 * t83 + t85 * t91;
t62 = atan2(t65, t74);
t59 = sin(t62);
t60 = cos(t62);
t53 = t59 * t65 + t60 * t74;
t52 = 0.1e1 / t53 ^ 2;
t72 = t79 * t89 + t81 * t84;
t63 = -t72 * t85 + t79 * t94;
t98 = t52 * t63;
t92 = t80 * t85;
t64 = t72 * t83 + t79 * t92;
t90 = t82 * t84;
t73 = -t79 * t90 + t81 * t86;
t78 = qJ(5) + qJ(6);
t76 = sin(t78);
t77 = cos(t78);
t58 = t64 * t77 + t73 * t76;
t56 = 0.1e1 / t58 ^ 2;
t57 = t64 * t76 - t73 * t77;
t97 = t56 * t57;
t69 = 0.1e1 / t74 ^ 2;
t96 = t65 * t69;
t95 = t73 * t83;
t93 = t80 * t84;
t88 = t56 * t57 ^ 2 + 0.1e1;
t87 = -t59 * t74 + t60 * t65;
t75 = t82 * t85 - t83 * t91;
t71 = t79 * t86 + t81 * t90;
t68 = 0.1e1 / t74;
t66 = -t70 * t83 + t81 * t92;
t61 = 0.1e1 / (t65 ^ 2 * t69 + 0.1e1);
t55 = 0.1e1 / t58;
t54 = 0.1e1 / t88;
t51 = 0.1e1 / t53;
t50 = 0.1e1 / (t52 * t63 ^ 2 + 0.1e1);
t49 = (t68 * t71 + t93 * t96) * t85 * t61;
t48 = (t66 * t68 - t75 * t96) * t61;
t47 = t88 * t54;
t1 = [0, t49, 0, t48, 0, 0; 0 (-t73 * t85 * t51 - ((t59 * t71 - t60 * t93) * t85 + t87 * t49) * t98) * t50, 0 (t64 * t51 - (t48 * t87 + t59 * t66 + t60 * t75) * t98) * t50, 0, 0; 0 ((t72 * t77 + t76 * t95) * t55 - (-t72 * t76 + t77 * t95) * t97) * t54, 0 (-t55 * t76 + t77 * t97) * t63 * t54, t47, t47;];
Ja_rot  = t1;
