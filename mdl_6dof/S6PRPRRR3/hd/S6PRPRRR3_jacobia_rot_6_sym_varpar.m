% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:48
% EndTime: 2019-02-26 19:54:48
% DurationCPUTime: 0.15s
% Computational Cost: add. (1407->31), mult. (1387->79), div. (95->9), fcn. (1976->13), ass. (0->54)
t93 = sin(pkin(6));
t94 = cos(pkin(11));
t106 = t93 * t94;
t95 = cos(pkin(6));
t97 = sin(qJ(2));
t103 = t95 * t97;
t92 = sin(pkin(11));
t99 = cos(qJ(2));
t85 = t94 * t103 + t92 * t99;
t91 = pkin(12) + qJ(4) + qJ(5);
t89 = sin(t91);
t90 = cos(t91);
t75 = t90 * t106 + t85 * t89;
t105 = t93 * t97;
t82 = t89 * t105 - t95 * t90;
t70 = atan2(-t75, t82);
t67 = sin(t70);
t68 = cos(t70);
t65 = -t67 * t75 + t68 * t82;
t64 = 0.1e1 / t65 ^ 2;
t107 = t92 * t93;
t87 = -t92 * t103 + t94 * t99;
t78 = -t90 * t107 + t87 * t89;
t112 = t64 * t78;
t102 = t95 * t99;
t86 = t92 * t102 + t94 * t97;
t96 = sin(qJ(6));
t109 = t86 * t96;
t79 = t89 * t107 + t87 * t90;
t98 = cos(qJ(6));
t74 = t79 * t98 + t109;
t72 = 0.1e1 / t74 ^ 2;
t108 = t86 * t98;
t73 = t79 * t96 - t108;
t111 = t72 * t73;
t81 = 0.1e1 / t82 ^ 2;
t110 = t75 * t81;
t104 = t93 * t99;
t101 = t73 ^ 2 * t72 + 0.1e1;
t100 = -t67 * t82 - t68 * t75;
t84 = t94 * t102 - t92 * t97;
t83 = t90 * t105 + t95 * t89;
t80 = 0.1e1 / t82;
t77 = -t89 * t106 + t85 * t90;
t71 = 0.1e1 / t74;
t69 = 0.1e1 / (t75 ^ 2 * t81 + 0.1e1);
t66 = 0.1e1 / t101;
t63 = 0.1e1 / t65;
t62 = 0.1e1 / (t78 ^ 2 * t64 + 0.1e1);
t61 = (t104 * t110 - t80 * t84) * t89 * t69;
t60 = (t83 * t110 - t77 * t80) * t69;
t59 = (t98 * t111 - t71 * t96) * t78 * t66;
t58 = (t79 * t63 - (t100 * t60 - t67 * t77 + t68 * t83) * t112) * t62;
t1 = [0, t61, 0, t60, t60, 0; 0 (-t86 * t89 * t63 - ((t68 * t104 - t67 * t84) * t89 + t100 * t61) * t112) * t62, 0, t58, t58, 0; 0 ((-t90 * t109 - t87 * t98) * t71 - (-t90 * t108 + t87 * t96) * t111) * t66, 0, t59, t59, t101 * t66;];
Ja_rot  = t1;
