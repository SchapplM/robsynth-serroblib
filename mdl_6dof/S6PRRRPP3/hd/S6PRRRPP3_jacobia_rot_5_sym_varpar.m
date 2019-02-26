% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPP3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:01
% EndTime: 2019-02-26 20:10:01
% DurationCPUTime: 0.22s
% Computational Cost: add. (707->41), mult. (1976->104), div. (87->9), fcn. (2786->13), ass. (0->56)
t87 = sin(pkin(6));
t91 = sin(qJ(3));
t103 = t87 * t91;
t89 = cos(pkin(6));
t92 = sin(qJ(2));
t101 = t89 * t92;
t86 = sin(pkin(10));
t88 = cos(pkin(10));
t95 = cos(qJ(2));
t80 = t88 * t101 + t86 * t95;
t94 = cos(qJ(3));
t71 = -t88 * t103 + t80 * t94;
t90 = sin(qJ(4));
t93 = cos(qJ(4));
t100 = t89 * t95;
t96 = t88 * t100 - t86 * t92;
t62 = t71 * t90 + t96 * t93;
t102 = t87 * t94;
t84 = t92 * t102 + t89 * t91;
t76 = t87 * t95 * t93 + t84 * t90;
t61 = atan2(-t62, t76);
t58 = sin(t61);
t59 = cos(t61);
t56 = -t58 * t62 + t59 * t76;
t55 = 0.1e1 / t56 ^ 2;
t81 = t86 * t100 + t88 * t92;
t104 = t81 * t93;
t82 = -t86 * t101 + t88 * t95;
t73 = t86 * t103 + t82 * t94;
t65 = t73 * t90 - t104;
t107 = t55 * t65;
t75 = 0.1e1 / t76 ^ 2;
t106 = t62 * t75;
t66 = t73 * t93 + t81 * t90;
t72 = -t86 * t102 + t82 * t91;
t69 = 0.1e1 / t72 ^ 2;
t105 = t66 * t69;
t99 = t90 * t94;
t98 = t90 * t95;
t97 = -t58 * t76 - t59 * t62;
t83 = -t92 * t103 + t89 * t94;
t78 = (-t92 * t93 + t94 * t98) * t87;
t77 = t84 * t93 - t87 * t98;
t74 = 0.1e1 / t76;
t70 = -t88 * t102 - t80 * t91;
t68 = 0.1e1 / t72;
t67 = -t80 * t93 + t96 * t99;
t64 = t71 * t93 - t96 * t90;
t60 = 0.1e1 / (t62 ^ 2 * t75 + 0.1e1);
t57 = 0.1e1 / (t66 ^ 2 * t69 + 0.1e1);
t54 = 0.1e1 / t56;
t53 = 0.1e1 / (t65 ^ 2 * t55 + 0.1e1);
t52 = (t83 * t106 - t70 * t74) * t90 * t60;
t51 = (t78 * t106 - t67 * t74) * t60;
t50 = (t77 * t106 - t64 * t74) * t60;
t1 = [0, t51, t52, t50, 0, 0; 0 ((-t81 * t99 - t82 * t93) * t54 - (t97 * t51 - t58 * t67 + t59 * t78) * t107) * t53 (-t72 * t90 * t54 - ((-t58 * t70 + t59 * t83) * t90 + t97 * t52) * t107) * t53 (t66 * t54 - (t97 * t50 - t58 * t64 + t59 * t77) * t107) * t53, 0, 0; 0 ((-t94 * t104 + t82 * t90) * t68 + t81 * t91 * t105) * t57 (-t68 * t72 * t93 - t73 * t105) * t57, -t65 * t68 * t57, 0, 0;];
Ja_rot  = t1;
