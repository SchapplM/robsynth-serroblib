% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR9_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:13
% EndTime: 2019-02-26 22:35:13
% DurationCPUTime: 0.24s
% Computational Cost: add. (1140->38), mult. (1655->89), div. (110->9), fcn. (2373->13), ass. (0->57)
t101 = cos(qJ(1));
t95 = sin(pkin(6));
t106 = t101 * t95;
t100 = cos(qJ(2));
t99 = sin(qJ(1));
t104 = t99 * t100;
t98 = sin(qJ(2));
t105 = t101 * t98;
t97 = cos(pkin(6));
t86 = t105 * t97 + t104;
t93 = qJ(3) + qJ(4);
t91 = sin(t93);
t92 = cos(t93);
t75 = t106 * t92 + t86 * t91;
t110 = t95 * t98;
t83 = t110 * t91 - t92 * t97;
t74 = atan2(-t75, t83);
t71 = sin(t74);
t72 = cos(t74);
t65 = -t71 * t75 + t72 * t83;
t64 = 0.1e1 / t65 ^ 2;
t109 = t95 * t99;
t103 = t101 * t100;
t108 = t99 * t98;
t88 = -t108 * t97 + t103;
t79 = -t109 * t92 + t88 * t91;
t117 = t64 * t79;
t87 = t104 * t97 + t105;
t94 = sin(pkin(12));
t112 = t87 * t94;
t80 = t109 * t91 + t88 * t92;
t96 = cos(pkin(12));
t70 = t80 * t96 + t112;
t68 = 0.1e1 / t70 ^ 2;
t111 = t87 * t96;
t69 = t80 * t94 - t111;
t116 = t68 * t69;
t115 = t72 * t75;
t82 = 0.1e1 / t83 ^ 2;
t114 = t75 * t82;
t113 = t79 ^ 2 * t64;
t107 = t100 * t95;
t77 = -t106 * t91 + t86 * t92;
t102 = -t71 * t83 - t115;
t85 = t103 * t97 - t108;
t84 = t110 * t92 + t91 * t97;
t81 = 0.1e1 / t83;
t73 = 0.1e1 / (t75 ^ 2 * t82 + 0.1e1);
t67 = 0.1e1 / t70;
t66 = 0.1e1 / (t68 * t69 ^ 2 + 0.1e1);
t63 = 0.1e1 / t65;
t62 = 0.1e1 / (0.1e1 + t113);
t61 = (t107 * t114 - t81 * t85) * t91 * t73;
t60 = (t114 * t84 - t77 * t81) * t73;
t59 = (t116 * t96 - t67 * t94) * t79 * t66;
t58 = (t80 * t63 - (t102 * t60 - t71 * t77 + t72 * t84) * t117) * t62;
t1 = [-t79 * t81 * t73, t61, t60, t60, 0, 0; (-t75 * t63 - (-t71 + (t115 * t81 + t71) * t73) * t113) * t62 (-t87 * t91 * t63 - ((t107 * t72 - t71 * t85) * t91 + t102 * t61) * t117) * t62, t58, t58, 0, 0; ((-t77 * t94 - t85 * t96) * t67 - (-t77 * t96 + t85 * t94) * t116) * t66 ((-t112 * t92 - t88 * t96) * t67 - (-t111 * t92 + t88 * t94) * t116) * t66, t59, t59, 0, 0;];
Ja_rot  = t1;
