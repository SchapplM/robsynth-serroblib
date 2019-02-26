% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR9_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:09
% EndTime: 2019-02-26 22:08:10
% DurationCPUTime: 0.27s
% Computational Cost: add. (666->39), mult. (1860->96), div. (80->9), fcn. (2626->13), ass. (0->58)
t96 = cos(qJ(2));
t97 = cos(qJ(1));
t102 = t97 * t96;
t93 = sin(qJ(2));
t94 = sin(qJ(1));
t105 = t94 * t93;
t91 = cos(pkin(6));
t101 = -t91 * t102 + t105;
t89 = sin(pkin(6));
t106 = t89 * t97;
t103 = t97 * t93;
t104 = t94 * t96;
t84 = t91 * t103 + t104;
t92 = sin(qJ(3));
t95 = cos(qJ(3));
t77 = t92 * t106 - t84 * t95;
t88 = sin(pkin(11));
t90 = cos(pkin(11));
t118 = t101 * t90 + t77 * t88;
t107 = t89 * t95;
t74 = (t93 * t107 + t91 * t92) * t88 + t89 * t96 * t90;
t63 = atan2(t118, t74);
t60 = sin(t63);
t61 = cos(t63);
t59 = t118 * t60 + t61 * t74;
t58 = 0.1e1 / t59 ^ 2;
t85 = t91 * t104 + t103;
t110 = t85 * t90;
t108 = t89 * t92;
t86 = -t91 * t105 + t102;
t79 = t94 * t108 + t86 * t95;
t69 = t79 * t88 - t110;
t117 = t58 * t69;
t116 = t61 * t118;
t70 = t79 * t90 + t85 * t88;
t66 = 0.1e1 / t70 ^ 2;
t78 = t94 * t107 - t86 * t92;
t115 = t66 * t78;
t73 = 0.1e1 / t74 ^ 2;
t114 = t118 * t73;
t113 = t69 ^ 2 * t58;
t112 = t78 ^ 2 * t66;
t109 = t88 * t95;
t100 = t101 * t88;
t99 = -t60 * t74 + t116;
t98 = t95 * t106 + t84 * t92;
t83 = -t93 * t108 + t91 * t95;
t80 = (t96 * t109 - t90 * t93) * t89;
t72 = 0.1e1 / t74;
t71 = -t95 * t100 - t84 * t90;
t65 = 0.1e1 / t70;
t64 = 0.1e1 / (0.1e1 + t112);
t62 = 0.1e1 / (t118 ^ 2 * t73 + 0.1e1);
t57 = 0.1e1 / t59;
t56 = 0.1e1 / (0.1e1 + t113);
t55 = (-t83 * t114 + t72 * t98) * t88 * t62;
t54 = (-t80 * t114 - t71 * t72) * t62;
t1 = [-t69 * t72 * t62, t54, t55, 0, 0, 0; (t118 * t57 - (-t60 + (-t72 * t116 + t60) * t62) * t113) * t56 ((-t85 * t109 - t86 * t90) * t57 - (t99 * t54 - t60 * t71 + t61 * t80) * t117) * t56 (t78 * t88 * t57 - ((t60 * t98 + t61 * t83) * t88 + t99 * t55) * t117) * t56, 0, 0, 0; (t98 * t65 - (t77 * t90 - t100) * t115) * t64 (t85 * t92 * t65 - (-t95 * t110 + t86 * t88) * t115) * t64 (-t90 * t112 - t65 * t79) * t64, 0, 0, 0;];
Ja_rot  = t1;
