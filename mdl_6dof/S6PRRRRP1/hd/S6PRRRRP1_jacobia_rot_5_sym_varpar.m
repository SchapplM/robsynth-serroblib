% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRP1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:12
% EndTime: 2019-02-26 20:15:13
% DurationCPUTime: 0.15s
% Computational Cost: add. (948->31), mult. (1387->79), div. (95->9), fcn. (1976->13), ass. (0->54)
t92 = sin(pkin(6));
t93 = cos(pkin(11));
t105 = t92 * t93;
t94 = cos(pkin(6));
t96 = sin(qJ(2));
t102 = t94 * t96;
t91 = sin(pkin(11));
t98 = cos(qJ(2));
t84 = t93 * t102 + t91 * t98;
t90 = qJ(3) + qJ(4);
t88 = sin(t90);
t89 = cos(t90);
t74 = t89 * t105 + t84 * t88;
t104 = t92 * t96;
t81 = t88 * t104 - t94 * t89;
t73 = atan2(-t74, t81);
t70 = sin(t73);
t71 = cos(t73);
t64 = -t70 * t74 + t71 * t81;
t63 = 0.1e1 / t64 ^ 2;
t106 = t91 * t92;
t86 = -t91 * t102 + t93 * t98;
t77 = -t89 * t106 + t86 * t88;
t111 = t63 * t77;
t101 = t94 * t98;
t85 = t91 * t101 + t93 * t96;
t95 = sin(qJ(5));
t108 = t85 * t95;
t78 = t88 * t106 + t86 * t89;
t97 = cos(qJ(5));
t69 = t78 * t97 + t108;
t67 = 0.1e1 / t69 ^ 2;
t107 = t85 * t97;
t68 = t78 * t95 - t107;
t110 = t67 * t68;
t80 = 0.1e1 / t81 ^ 2;
t109 = t74 * t80;
t103 = t92 * t98;
t100 = t68 ^ 2 * t67 + 0.1e1;
t99 = -t70 * t81 - t71 * t74;
t83 = t93 * t101 - t91 * t96;
t82 = t89 * t104 + t94 * t88;
t79 = 0.1e1 / t81;
t76 = -t88 * t105 + t84 * t89;
t72 = 0.1e1 / (t74 ^ 2 * t80 + 0.1e1);
t66 = 0.1e1 / t69;
t65 = 0.1e1 / t100;
t62 = 0.1e1 / t64;
t61 = 0.1e1 / (t77 ^ 2 * t63 + 0.1e1);
t60 = (t103 * t109 - t79 * t83) * t88 * t72;
t59 = (t82 * t109 - t76 * t79) * t72;
t58 = (t97 * t110 - t66 * t95) * t77 * t65;
t57 = (t78 * t62 - (t99 * t59 - t70 * t76 + t71 * t82) * t111) * t61;
t1 = [0, t60, t59, t59, 0, 0; 0 (-t85 * t88 * t62 - ((t71 * t103 - t70 * t83) * t88 + t99 * t60) * t111) * t61, t57, t57, 0, 0; 0 ((-t89 * t108 - t86 * t97) * t66 - (-t89 * t107 + t86 * t95) * t110) * t65, t58, t58, t100 * t65, 0;];
Ja_rot  = t1;
