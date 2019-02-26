% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:50
% EndTime: 2019-02-26 22:05:50
% DurationCPUTime: 0.22s
% Computational Cost: add. (944->39), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->56)
t88 = sin(pkin(6));
t93 = cos(qJ(1));
t100 = t88 * t93;
t89 = cos(pkin(6));
t90 = sin(qJ(2));
t97 = t93 * t90;
t91 = sin(qJ(1));
t92 = cos(qJ(2));
t98 = t91 * t92;
t77 = t89 * t97 + t98;
t87 = qJ(3) + pkin(11);
t83 = sin(t87);
t85 = cos(t87);
t66 = t85 * t100 + t77 * t83;
t103 = t88 * t90;
t74 = t83 * t103 - t89 * t85;
t65 = atan2(-t66, t74);
t62 = sin(t65);
t63 = cos(t65);
t56 = -t62 * t66 + t63 * t74;
t55 = 0.1e1 / t56 ^ 2;
t102 = t88 * t91;
t96 = t93 * t92;
t99 = t91 * t90;
t79 = -t89 * t99 + t96;
t70 = -t85 * t102 + t79 * t83;
t109 = t55 * t70;
t71 = t83 * t102 + t79 * t85;
t78 = t89 * t98 + t97;
t86 = pkin(12) + qJ(6);
t82 = sin(t86);
t84 = cos(t86);
t61 = t71 * t84 + t78 * t82;
t59 = 0.1e1 / t61 ^ 2;
t60 = t71 * t82 - t78 * t84;
t108 = t59 * t60;
t107 = t63 * t66;
t73 = 0.1e1 / t74 ^ 2;
t106 = t66 * t73;
t105 = t70 ^ 2 * t55;
t104 = t78 * t85;
t101 = t88 * t92;
t95 = t60 ^ 2 * t59 + 0.1e1;
t68 = -t83 * t100 + t77 * t85;
t94 = -t62 * t74 - t107;
t76 = t89 * t96 - t99;
t75 = t85 * t103 + t89 * t83;
t72 = 0.1e1 / t74;
t64 = 0.1e1 / (t66 ^ 2 * t73 + 0.1e1);
t58 = 0.1e1 / t61;
t57 = 0.1e1 / t95;
t54 = 0.1e1 / t56;
t53 = 0.1e1 / (0.1e1 + t105);
t52 = (t101 * t106 - t72 * t76) * t83 * t64;
t51 = (t75 * t106 - t68 * t72) * t64;
t1 = [-t70 * t72 * t64, t52, t51, 0, 0, 0; (-t66 * t54 - (-t62 + (t72 * t107 + t62) * t64) * t105) * t53 (-t78 * t83 * t54 - ((t63 * t101 - t62 * t76) * t83 + t94 * t52) * t109) * t53 (t71 * t54 - (t94 * t51 - t62 * t68 + t63 * t75) * t109) * t53, 0, 0, 0; ((-t68 * t82 - t76 * t84) * t58 - (-t68 * t84 + t76 * t82) * t108) * t57 ((-t82 * t104 - t79 * t84) * t58 - (-t84 * t104 + t79 * t82) * t108) * t57 (t84 * t108 - t82 * t58) * t70 * t57, 0, 0, t95 * t57;];
Ja_rot  = t1;
