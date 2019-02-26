% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPP1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:44
% EndTime: 2019-02-26 20:08:44
% DurationCPUTime: 0.23s
% Computational Cost: add. (1132->41), mult. (1993->101), div. (87->9), fcn. (2807->13), ass. (0->58)
t101 = cos(pkin(10));
t104 = sin(qJ(2));
t102 = cos(pkin(6));
t106 = cos(qJ(2));
t109 = t102 * t106;
t99 = sin(pkin(10));
t107 = t101 * t109 - t99 * t104;
t105 = cos(qJ(3));
t100 = sin(pkin(6));
t103 = sin(qJ(3));
t113 = t100 * t103;
t110 = t102 * t104;
t91 = t101 * t110 + t99 * t106;
t86 = -t101 * t113 + t91 * t105;
t98 = qJ(4) + pkin(11);
t96 = sin(t98);
t97 = cos(t98);
t74 = t107 * t97 + t86 * t96;
t111 = t100 * t106;
t112 = t100 * t105;
t95 = t102 * t103 + t104 * t112;
t82 = t97 * t111 + t95 * t96;
t70 = atan2(-t74, t82);
t67 = sin(t70);
t68 = cos(t70);
t66 = -t67 * t74 + t68 * t82;
t65 = 0.1e1 / t66 ^ 2;
t92 = t101 * t104 + t99 * t109;
t115 = t92 * t97;
t93 = t101 * t106 - t99 * t110;
t88 = t93 * t105 + t99 * t113;
t77 = t88 * t96 - t115;
t119 = t65 * t77;
t78 = t88 * t97 + t92 * t96;
t73 = 0.1e1 / t78 ^ 2;
t87 = -t93 * t103 + t99 * t112;
t118 = t73 * t87;
t81 = 0.1e1 / t82 ^ 2;
t117 = t74 * t81;
t116 = t87 ^ 2 * t73;
t114 = t105 * t96;
t108 = -t67 * t82 - t68 * t74;
t94 = t102 * t105 - t104 * t113;
t89 = (-t104 * t97 + t106 * t114) * t100;
t85 = -t101 * t112 - t91 * t103;
t83 = -t96 * t111 + t95 * t97;
t80 = 0.1e1 / t82;
t79 = t107 * t114 - t91 * t97;
t76 = -t107 * t96 + t86 * t97;
t72 = 0.1e1 / t78;
t71 = 0.1e1 / (0.1e1 + t116);
t69 = 0.1e1 / (t74 ^ 2 * t81 + 0.1e1);
t64 = 0.1e1 / t66;
t63 = 0.1e1 / (t77 ^ 2 * t65 + 0.1e1);
t62 = (t94 * t117 - t80 * t85) * t96 * t69;
t61 = (t89 * t117 - t79 * t80) * t69;
t60 = (t83 * t117 - t76 * t80) * t69;
t1 = [0, t61, t62, t60, 0, 0; 0 ((-t92 * t114 - t93 * t97) * t64 - (t108 * t61 - t67 * t79 + t68 * t89) * t119) * t63 (t87 * t96 * t64 - ((-t67 * t85 + t68 * t94) * t96 + t108 * t62) * t119) * t63 (t78 * t64 - (t108 * t60 - t67 * t76 + t68 * t83) * t119) * t63, 0, 0; 0 (t92 * t103 * t72 - (-t105 * t115 + t93 * t96) * t118) * t71 (-t97 * t116 - t72 * t88) * t71, t77 * t71 * t118, 0, 0;];
Ja_rot  = t1;
