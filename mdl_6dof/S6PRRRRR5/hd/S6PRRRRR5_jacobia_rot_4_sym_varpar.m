% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRR5_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_jacobia_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:20
% EndTime: 2019-02-26 20:21:20
% DurationCPUTime: 0.22s
% Computational Cost: add. (607->40), mult. (1825->97), div. (65->9), fcn. (2498->15), ass. (0->60)
t87 = sin(pkin(7));
t88 = sin(pkin(6));
t109 = t88 * t87;
t89 = cos(pkin(13));
t90 = cos(pkin(7));
t91 = cos(pkin(6));
t97 = cos(qJ(2));
t106 = t91 * t97;
t86 = sin(pkin(13));
t94 = sin(qJ(2));
t99 = t89 * t106 - t86 * t94;
t116 = -t89 * t109 + t99 * t90;
t107 = t91 * t94;
t81 = t89 * t107 + t86 * t97;
t93 = sin(qJ(3));
t96 = cos(qJ(3));
t66 = -t116 * t96 + t81 * t93;
t108 = t90 * t96;
t110 = t87 * t91;
t75 = -t96 * t110 + (-t97 * t108 + t93 * t94) * t88;
t65 = atan2(-t66, t75);
t62 = sin(t65);
t63 = cos(t65);
t56 = -t62 * t66 + t63 * t75;
t55 = 0.1e1 / t56 ^ 2;
t103 = t86 * t109;
t83 = -t86 * t107 + t89 * t97;
t111 = t83 * t93;
t82 = -t86 * t106 - t89 * t94;
t69 = -t96 * t103 - t82 * t108 + t111;
t115 = t55 * t69;
t70 = t83 * t96 + (t82 * t90 + t103) * t93;
t77 = t86 * t88 * t90 - t82 * t87;
t92 = sin(qJ(4));
t95 = cos(qJ(4));
t61 = t70 * t95 + t77 * t92;
t59 = 0.1e1 / t61 ^ 2;
t60 = t70 * t92 - t77 * t95;
t114 = t59 * t60;
t74 = 0.1e1 / t75 ^ 2;
t113 = t66 * t74;
t112 = t83 * t87;
t105 = t93 * t97;
t104 = t94 * t96;
t101 = t60 ^ 2 * t59 + 0.1e1;
t100 = -t62 * t75 - t63 * t66;
t80 = (t90 * t104 + t105) * t88;
t76 = t93 * t110 + (t90 * t105 + t104) * t88;
t73 = 0.1e1 / t75;
t72 = -t90 * t111 + t82 * t96;
t71 = t81 * t108 + t99 * t93;
t68 = t116 * t93 + t81 * t96;
t64 = 0.1e1 / (t66 ^ 2 * t74 + 0.1e1);
t58 = 0.1e1 / t61;
t57 = 0.1e1 / t101;
t54 = 0.1e1 / t56;
t53 = 0.1e1 / (t69 ^ 2 * t55 + 0.1e1);
t52 = (t80 * t113 - t71 * t73) * t64;
t51 = (t76 * t113 - t68 * t73) * t64;
t1 = [0, t52, t51, 0, 0, 0; 0 ((t83 * t108 + t82 * t93) * t54 - (t100 * t52 - t62 * t71 + t63 * t80) * t115) * t53 (t70 * t54 - (t100 * t51 - t62 * t68 + t63 * t76) * t115) * t53, 0, 0, 0; 0 ((-t95 * t112 + t72 * t92) * t58 - (t92 * t112 + t72 * t95) * t114) * t57 (t95 * t114 - t92 * t58) * t69 * t57, t101 * t57, 0, 0;];
Ja_rot  = t1;
