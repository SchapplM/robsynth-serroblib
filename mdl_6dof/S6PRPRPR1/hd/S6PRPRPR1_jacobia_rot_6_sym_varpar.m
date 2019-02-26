% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:46:27
% EndTime: 2019-02-26 19:46:28
% DurationCPUTime: 0.18s
% Computational Cost: add. (939->33), mult. (1727->83), div. (65->9), fcn. (2425->15), ass. (0->56)
t93 = sin(pkin(6));
t95 = cos(pkin(10));
t107 = t93 * t95;
t100 = cos(qJ(2));
t91 = sin(pkin(11));
t94 = cos(pkin(11));
t98 = sin(qJ(2));
t103 = t100 * t94 - t98 * t91;
t102 = t100 * t91 + t98 * t94;
t96 = cos(pkin(6));
t83 = t102 * t96;
t92 = sin(pkin(10));
t72 = t103 * t92 + t95 * t83;
t90 = qJ(4) + pkin(12);
t88 = sin(t90);
t89 = cos(t90);
t66 = t89 * t107 + t72 * t88;
t82 = t102 * t93;
t78 = t82 * t88 - t96 * t89;
t65 = atan2(-t66, t78);
t62 = sin(t65);
t63 = cos(t65);
t56 = -t62 * t66 + t63 * t78;
t55 = 0.1e1 / t56 ^ 2;
t104 = t103 * t95 - t92 * t83;
t108 = t92 * t93;
t69 = t104 * t88 - t89 * t108;
t113 = t55 * t69;
t101 = t103 * t96;
t74 = -t92 * t101 - t102 * t95;
t97 = sin(qJ(6));
t110 = t74 * t97;
t70 = t104 * t89 + t88 * t108;
t99 = cos(qJ(6));
t61 = t70 * t99 - t110;
t59 = 0.1e1 / t61 ^ 2;
t109 = t74 * t99;
t60 = t70 * t97 + t109;
t112 = t59 * t60;
t77 = 0.1e1 / t78 ^ 2;
t111 = t66 * t77;
t106 = t60 ^ 2 * t59 + 0.1e1;
t105 = -t62 * t78 - t63 * t66;
t81 = t103 * t93;
t79 = t82 * t89 + t96 * t88;
t76 = 0.1e1 / t78;
t71 = t95 * t101 - t102 * t92;
t68 = -t88 * t107 + t72 * t89;
t64 = 0.1e1 / (t66 ^ 2 * t77 + 0.1e1);
t58 = 0.1e1 / t61;
t57 = 0.1e1 / t106;
t54 = 0.1e1 / t56;
t53 = 0.1e1 / (t69 ^ 2 * t55 + 0.1e1);
t52 = (t81 * t111 - t71 * t76) * t88 * t64;
t51 = (t79 * t111 - t68 * t76) * t64;
t1 = [0, t52, 0, t51, 0, 0; 0 (t74 * t88 * t54 - ((-t62 * t71 + t63 * t81) * t88 + t105 * t52) * t113) * t53, 0 (t70 * t54 - (t105 * t51 - t62 * t68 + t63 * t79) * t113) * t53, 0, 0; 0 ((-t104 * t99 + t89 * t110) * t58 - (t104 * t97 + t89 * t109) * t112) * t57, 0 (t99 * t112 - t58 * t97) * t69 * t57, 0, t106 * t57;];
Ja_rot  = t1;
