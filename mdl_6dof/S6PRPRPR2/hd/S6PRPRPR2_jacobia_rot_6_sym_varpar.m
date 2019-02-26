% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR2
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
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:03
% EndTime: 2019-02-26 19:47:03
% DurationCPUTime: 0.21s
% Computational Cost: add. (680->33), mult. (1727->84), div. (65->9), fcn. (2425->15), ass. (0->55)
t100 = cos(qJ(4));
t94 = sin(pkin(6));
t108 = t100 * t94;
t101 = cos(qJ(2));
t92 = sin(pkin(11));
t95 = cos(pkin(11));
t99 = sin(qJ(2));
t104 = t101 * t95 - t99 * t92;
t103 = t101 * t92 + t99 * t95;
t97 = cos(pkin(6));
t84 = t103 * t97;
t93 = sin(pkin(10));
t96 = cos(pkin(10));
t73 = t104 * t93 + t96 * t84;
t98 = sin(qJ(4));
t67 = t96 * t108 + t73 * t98;
t83 = t103 * t94;
t79 = -t97 * t100 + t83 * t98;
t66 = atan2(-t67, t79);
t63 = sin(t66);
t64 = cos(t66);
t57 = -t63 * t67 + t64 * t79;
t56 = 0.1e1 / t57 ^ 2;
t105 = t104 * t96 - t93 * t84;
t70 = t105 * t98 - t93 * t108;
t113 = t56 * t70;
t110 = t94 * t98;
t71 = t100 * t105 + t93 * t110;
t102 = t104 * t97;
t75 = -t93 * t102 - t103 * t96;
t91 = pkin(12) + qJ(6);
t89 = sin(t91);
t90 = cos(t91);
t62 = t71 * t90 - t75 * t89;
t60 = 0.1e1 / t62 ^ 2;
t61 = t71 * t89 + t75 * t90;
t112 = t60 * t61;
t78 = 0.1e1 / t79 ^ 2;
t111 = t67 * t78;
t109 = t100 * t75;
t107 = t61 ^ 2 * t60 + 0.1e1;
t106 = -t63 * t79 - t64 * t67;
t82 = t104 * t94;
t80 = t83 * t100 + t97 * t98;
t77 = 0.1e1 / t79;
t72 = t96 * t102 - t103 * t93;
t69 = t73 * t100 - t96 * t110;
t65 = 0.1e1 / (t67 ^ 2 * t78 + 0.1e1);
t59 = 0.1e1 / t62;
t58 = 0.1e1 / t107;
t55 = 0.1e1 / t57;
t54 = 0.1e1 / (t70 ^ 2 * t56 + 0.1e1);
t53 = (t82 * t111 - t72 * t77) * t98 * t65;
t52 = (t80 * t111 - t69 * t77) * t65;
t1 = [0, t53, 0, t52, 0, 0; 0 (t75 * t98 * t55 - ((-t63 * t72 + t64 * t82) * t98 + t106 * t53) * t113) * t54, 0 (t71 * t55 - (t106 * t52 - t63 * t69 + t64 * t80) * t113) * t54, 0, 0; 0 ((-t105 * t90 + t89 * t109) * t59 - (t105 * t89 + t90 * t109) * t112) * t58, 0 (t90 * t112 - t59 * t89) * t70 * t58, 0, t107 * t58;];
Ja_rot  = t1;
