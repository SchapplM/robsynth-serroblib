% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6PRRRPR6_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:13:25
% EndTime: 2019-02-26 20:13:25
% DurationCPUTime: 0.22s
% Computational Cost: add. (714->40), mult. (1993->102), div. (87->9), fcn. (2807->13), ass. (0->57)
t89 = sin(pkin(6));
t93 = sin(qJ(3));
t105 = t89 * t93;
t91 = cos(pkin(6));
t94 = sin(qJ(2));
t103 = t91 * t94;
t88 = sin(pkin(11));
t90 = cos(pkin(11));
t97 = cos(qJ(2));
t83 = t90 * t103 + t88 * t97;
t96 = cos(qJ(3));
t74 = -t90 * t105 + t83 * t96;
t92 = sin(qJ(4));
t95 = cos(qJ(4));
t102 = t91 * t97;
t98 = t90 * t102 - t88 * t94;
t66 = t74 * t92 + t98 * t95;
t104 = t89 * t96;
t87 = t94 * t104 + t91 * t93;
t79 = t89 * t97 * t95 + t87 * t92;
t63 = atan2(-t66, t79);
t60 = sin(t63);
t61 = cos(t63);
t58 = -t60 * t66 + t61 * t79;
t57 = 0.1e1 / t58 ^ 2;
t84 = t88 * t102 + t90 * t94;
t106 = t84 * t95;
t85 = -t88 * t103 + t90 * t97;
t76 = t88 * t105 + t85 * t96;
t69 = t76 * t92 - t106;
t110 = t57 * t69;
t70 = t76 * t95 + t84 * t92;
t65 = 0.1e1 / t70 ^ 2;
t75 = t88 * t104 - t85 * t93;
t109 = t65 * t75;
t78 = 0.1e1 / t79 ^ 2;
t108 = t66 * t78;
t107 = t75 ^ 2 * t65;
t101 = t92 * t96;
t100 = t92 * t97;
t99 = -t60 * t79 - t61 * t66;
t86 = -t94 * t105 + t91 * t96;
t81 = (t96 * t100 - t94 * t95) * t89;
t80 = -t89 * t100 + t87 * t95;
t77 = 0.1e1 / t79;
t73 = -t90 * t104 - t83 * t93;
t71 = t98 * t101 - t83 * t95;
t68 = t74 * t95 - t98 * t92;
t64 = 0.1e1 / t70;
t62 = 0.1e1 / (t66 ^ 2 * t78 + 0.1e1);
t59 = 0.1e1 / (0.1e1 + t107);
t56 = 0.1e1 / t58;
t55 = 0.1e1 / (t69 ^ 2 * t57 + 0.1e1);
t54 = (t86 * t108 - t73 * t77) * t92 * t62;
t53 = (t81 * t108 - t71 * t77) * t62;
t52 = (t80 * t108 - t68 * t77) * t62;
t1 = [0, t53, t54, t52, 0, 0; 0 ((-t84 * t101 - t85 * t95) * t56 - (t99 * t53 - t60 * t71 + t61 * t81) * t110) * t55 (t75 * t92 * t56 - ((-t60 * t73 + t61 * t86) * t92 + t99 * t54) * t110) * t55 (t70 * t56 - (t99 * t52 - t60 * t68 + t61 * t80) * t110) * t55, 0, 0; 0 (t84 * t93 * t64 - (-t96 * t106 + t85 * t92) * t109) * t59 (-t95 * t107 - t64 * t76) * t59, t69 * t59 * t109, 0, 0;];
Ja_rot  = t1;
