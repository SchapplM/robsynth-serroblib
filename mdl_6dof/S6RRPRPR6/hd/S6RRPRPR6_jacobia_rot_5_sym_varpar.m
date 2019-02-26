% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR6_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:37
% EndTime: 2019-02-26 21:40:37
% DurationCPUTime: 0.22s
% Computational Cost: add. (675->32), mult. (1865->82), div. (77->9), fcn. (2642->13), ass. (0->50)
t85 = cos(pkin(6));
t82 = sin(pkin(11));
t84 = cos(pkin(11));
t87 = sin(qJ(2));
t90 = cos(qJ(2));
t93 = t90 * t82 + t87 * t84;
t76 = t93 * t85;
t77 = t87 * t82 - t90 * t84;
t88 = sin(qJ(1));
t91 = cos(qJ(1));
t65 = t91 * t76 - t88 * t77;
t86 = sin(qJ(4));
t89 = cos(qJ(4));
t83 = sin(pkin(6));
t95 = t83 * t91;
t57 = t65 * t86 + t89 * t95;
t75 = t93 * t83;
t71 = t75 * t86 - t85 * t89;
t56 = atan2(-t57, t71);
t53 = sin(t56);
t54 = cos(t56);
t51 = -t53 * t57 + t54 * t71;
t50 = 0.1e1 / t51 ^ 2;
t68 = -t88 * t76 - t91 * t77;
t96 = t83 * t88;
t60 = t68 * t86 - t89 * t96;
t101 = t50 * t60;
t100 = t54 * t57;
t70 = 0.1e1 / t71 ^ 2;
t99 = t57 * t70;
t98 = t60 ^ 2 * t50;
t61 = t68 * t89 + t86 * t96;
t92 = t77 * t85;
t66 = t88 * t92 - t91 * t93;
t63 = 0.1e1 / t66 ^ 2;
t97 = t61 * t63;
t59 = t65 * t89 - t86 * t95;
t94 = -t53 * t71 - t100;
t74 = t77 * t83;
t72 = t75 * t89 + t85 * t86;
t69 = 0.1e1 / t71;
t64 = -t88 * t93 - t91 * t92;
t62 = 0.1e1 / t66;
t55 = 0.1e1 / (t57 ^ 2 * t70 + 0.1e1);
t52 = 0.1e1 / (t61 ^ 2 * t63 + 0.1e1);
t49 = 0.1e1 / t51;
t48 = 0.1e1 / (0.1e1 + t98);
t47 = (-t64 * t69 - t74 * t99) * t86 * t55;
t46 = (-t59 * t69 + t72 * t99) * t55;
t1 = [-t60 * t69 * t55, t47, 0, t46, 0, 0; (-t57 * t49 - (-t53 + (t69 * t100 + t53) * t55) * t98) * t48 (t66 * t86 * t49 - ((-t53 * t64 - t54 * t74) * t86 + t94 * t47) * t101) * t48, 0 (t61 * t49 - (t94 * t46 - t53 * t59 + t54 * t72) * t101) * t48, 0, 0; (t59 * t62 - t64 * t97) * t52 (-t66 * t89 * t62 - t68 * t97) * t52, 0, t60 * t62 * t52, 0, 0;];
Ja_rot  = t1;
