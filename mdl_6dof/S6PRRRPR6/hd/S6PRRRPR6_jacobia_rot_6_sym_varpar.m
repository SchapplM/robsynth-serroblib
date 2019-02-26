% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:13:35
% EndTime: 2019-02-26 20:13:35
% DurationCPUTime: 0.19s
% Computational Cost: add. (537->36), mult. (1484->92), div. (71->9), fcn. (2087->15), ass. (0->56)
t88 = sin(pkin(6));
t97 = cos(qJ(3));
t104 = t88 * t97;
t90 = cos(pkin(6));
t94 = sin(qJ(2));
t102 = t90 * t94;
t87 = sin(pkin(11));
t89 = cos(pkin(11));
t98 = cos(qJ(2));
t82 = t89 * t102 + t87 * t98;
t93 = sin(qJ(3));
t74 = t89 * t104 + t82 * t93;
t105 = t88 * t93;
t85 = -t94 * t105 + t90 * t97;
t71 = atan2(t74, t85);
t68 = sin(t71);
t69 = cos(t71);
t61 = t68 * t74 + t69 * t85;
t60 = 0.1e1 / t61 ^ 2;
t84 = -t87 * t102 + t89 * t98;
t76 = t87 * t104 - t84 * t93;
t110 = t60 * t76;
t77 = t87 * t105 + t84 * t97;
t101 = t90 * t98;
t83 = t87 * t101 + t89 * t94;
t92 = sin(qJ(4));
t96 = cos(qJ(4));
t66 = t77 * t92 - t83 * t96;
t67 = t77 * t96 + t83 * t92;
t91 = sin(qJ(6));
t95 = cos(qJ(6));
t65 = t66 * t91 + t67 * t95;
t63 = 0.1e1 / t65 ^ 2;
t64 = -t66 * t95 + t67 * t91;
t109 = t63 * t64;
t108 = t64 ^ 2 * t63;
t80 = 0.1e1 / t85 ^ 2;
t107 = t74 * t80;
t106 = t83 * t97;
t103 = t88 * t98;
t100 = 0.1e1 + t108;
t99 = -t68 * t85 + t69 * t74;
t86 = -t94 * t104 - t90 * t93;
t81 = t89 * t101 - t87 * t94;
t79 = 0.1e1 / t85;
t75 = -t89 * t105 + t82 * t97;
t73 = -t96 * t106 + t84 * t92;
t72 = -t92 * t106 - t84 * t96;
t70 = 0.1e1 / (t74 ^ 2 * t80 + 0.1e1);
t62 = 0.1e1 / t65;
t59 = 0.1e1 / t61;
t58 = 0.1e1 / (t76 ^ 2 * t60 + 0.1e1);
t57 = (t103 * t107 + t79 * t81) * t93 * t70;
t56 = (-t86 * t107 + t75 * t79) * t70;
t55 = 0.1e1 / t100;
t1 = [0, t57, t56, 0, 0, 0; 0 (t83 * t93 * t59 - ((-t69 * t103 + t68 * t81) * t93 + t99 * t57) * t110) * t58 (-t77 * t59 - (t99 * t56 + t68 * t75 + t69 * t86) * t110) * t58, 0, 0, 0; 0 ((-t72 * t95 + t73 * t91) * t62 - (t72 * t91 + t73 * t95) * t109) * t55 ((t91 * t96 - t92 * t95) * t62 - (t91 * t92 + t95 * t96) * t109) * t55 * t76 (-t65 * t62 - t108) * t55, 0, t100 * t55;];
Ja_rot  = t1;
