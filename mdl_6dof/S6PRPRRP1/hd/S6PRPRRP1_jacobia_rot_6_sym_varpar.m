% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRP1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:50:18
% EndTime: 2019-02-26 19:50:18
% DurationCPUTime: 0.21s
% Computational Cost: add. (632->32), mult. (1727->83), div. (65->9), fcn. (2425->15), ass. (0->55)
t100 = sin(qJ(4));
t103 = cos(qJ(4));
t95 = sin(pkin(6));
t111 = t103 * t95;
t101 = sin(qJ(2));
t104 = cos(qJ(2));
t93 = sin(pkin(11));
t96 = cos(pkin(11));
t106 = t101 * t96 + t104 * t93;
t98 = cos(pkin(6));
t88 = t106 * t98;
t89 = t101 * t93 - t104 * t96;
t94 = sin(pkin(10));
t97 = cos(pkin(10));
t77 = t97 * t88 - t94 * t89;
t71 = t77 * t100 + t111 * t97;
t87 = t106 * t95;
t83 = t87 * t100 - t98 * t103;
t70 = atan2(-t71, t83);
t67 = sin(t70);
t68 = cos(t70);
t61 = -t67 * t71 + t68 * t83;
t60 = 0.1e1 / t61 ^ 2;
t107 = -t94 * t88 - t97 * t89;
t74 = t100 * t107 - t111 * t94;
t116 = t60 * t74;
t102 = cos(qJ(5));
t105 = t89 * t98;
t79 = t105 * t94 - t106 * t97;
t99 = sin(qJ(5));
t113 = t79 * t99;
t112 = t100 * t95;
t75 = t103 * t107 + t112 * t94;
t66 = t75 * t102 - t113;
t64 = 0.1e1 / t66 ^ 2;
t110 = t79 * t102;
t65 = t75 * t99 + t110;
t115 = t64 * t65;
t82 = 0.1e1 / t83 ^ 2;
t114 = t71 * t82;
t109 = t65 ^ 2 * t64 + 0.1e1;
t108 = -t67 * t83 - t68 * t71;
t86 = t89 * t95;
t84 = t98 * t100 + t87 * t103;
t81 = 0.1e1 / t83;
t76 = -t105 * t97 - t106 * t94;
t73 = t77 * t103 - t112 * t97;
t69 = 0.1e1 / (t71 ^ 2 * t82 + 0.1e1);
t63 = 0.1e1 / t66;
t62 = 0.1e1 / t109;
t59 = 0.1e1 / t61;
t58 = 0.1e1 / (t74 ^ 2 * t60 + 0.1e1);
t57 = (-t114 * t86 - t76 * t81) * t69 * t100;
t56 = (t114 * t84 - t73 * t81) * t69;
t1 = [0, t57, 0, t56, 0, 0; 0 (t79 * t100 * t59 - (t108 * t57 + (-t67 * t76 - t68 * t86) * t100) * t116) * t58, 0 (t75 * t59 - (t108 * t56 - t67 * t73 + t68 * t84) * t116) * t58, 0, 0; 0 ((-t102 * t107 + t103 * t113) * t63 - (t103 * t110 + t107 * t99) * t115) * t62, 0 (t102 * t115 - t63 * t99) * t74 * t62, t109 * t62, 0;];
Ja_rot  = t1;
