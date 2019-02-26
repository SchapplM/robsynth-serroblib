% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR10
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR10_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobia_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:49
% EndTime: 2019-02-26 21:19:50
% DurationCPUTime: 0.21s
% Computational Cost: add. (559->39), mult. (1669->87), div. (55->9), fcn. (2293->15), ass. (0->57)
t115 = sin(qJ(1));
t88 = sin(pkin(6));
t102 = t88 * t115;
t87 = sin(pkin(7));
t90 = cos(pkin(7));
t89 = cos(pkin(13));
t100 = t115 * t89;
t86 = sin(pkin(13));
t96 = cos(qJ(1));
t106 = t96 * t86;
t91 = cos(pkin(6));
t98 = t91 * t100 + t106;
t116 = -t87 * t102 + t98 * t90;
t101 = t115 * t86;
t105 = t96 * t89;
t83 = -t91 * t101 + t105;
t93 = sin(qJ(3));
t95 = cos(qJ(3));
t72 = -t116 * t93 + t83 * t95;
t78 = t90 * t102 + t98 * t87;
t92 = sin(qJ(4));
t94 = cos(qJ(4));
t62 = t72 * t94 + t78 * t92;
t60 = 0.1e1 / t62 ^ 2;
t61 = t72 * t92 - t78 * t94;
t114 = t60 * t61;
t109 = t88 * t96;
t104 = t87 * t109;
t107 = t90 * t95;
t82 = t91 * t106 + t100;
t111 = t82 * t93;
t81 = -t91 * t105 + t101;
t67 = t95 * t104 + t81 * t107 + t111;
t110 = t87 * t91;
t75 = -t95 * t110 + (-t89 * t107 + t86 * t93) * t88;
t66 = atan2(-t67, t75);
t64 = cos(t66);
t113 = t64 * t67;
t63 = sin(t66);
t57 = -t63 * t67 + t64 * t75;
t56 = 0.1e1 / t57 ^ 2;
t71 = t116 * t95 + t83 * t93;
t112 = t71 ^ 2 * t56;
t108 = t90 * t93;
t103 = t61 ^ 2 * t60 + 0.1e1;
t70 = t93 * t104 + t81 * t108 - t82 * t95;
t77 = t90 * t109 - t81 * t87;
t76 = t93 * t110 + (t89 * t108 + t86 * t95) * t88;
t74 = 0.1e1 / t75 ^ 2;
t73 = 0.1e1 / t75;
t65 = 0.1e1 / (t67 ^ 2 * t74 + 0.1e1);
t59 = 0.1e1 / t62;
t58 = 0.1e1 / t103;
t55 = 0.1e1 / t57;
t54 = 0.1e1 / (0.1e1 + t112);
t53 = (t67 * t74 * t76 + t70 * t73) * t65;
t1 = [-t71 * t73 * t65, 0, t53, 0, 0, 0; ((-t111 + (-t81 * t90 - t104) * t95) * t55 - (-t63 + (t73 * t113 + t63) * t65) * t112) * t54, 0 (t72 * t55 - (t63 * t70 + t64 * t76 + (-t63 * t75 - t113) * t53) * t71 * t56) * t54, 0, 0, 0; ((t70 * t92 - t77 * t94) * t59 - (t70 * t94 + t77 * t92) * t114) * t58, 0 (t94 * t114 - t92 * t59) * t71 * t58, t103 * t58, 0, 0;];
Ja_rot  = t1;
