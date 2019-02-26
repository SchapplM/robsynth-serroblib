% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:07
% EndTime: 2019-02-26 20:11:07
% DurationCPUTime: 0.14s
% Computational Cost: add. (1012->32), mult. (1387->80), div. (95->9), fcn. (1976->13), ass. (0->54)
t102 = sin(pkin(6));
t103 = cos(pkin(11));
t113 = t102 * t103;
t101 = sin(pkin(11));
t106 = cos(qJ(2));
t104 = cos(pkin(6));
t105 = sin(qJ(2));
t110 = t104 * t105;
t91 = t101 * t106 + t103 * t110;
t100 = qJ(3) + qJ(4);
t97 = sin(t100);
t98 = cos(t100);
t81 = t98 * t113 + t91 * t97;
t112 = t102 * t105;
t88 = -t104 * t98 + t97 * t112;
t80 = atan2(-t81, t88);
t77 = sin(t80);
t78 = cos(t80);
t71 = -t77 * t81 + t78 * t88;
t70 = 0.1e1 / t71 ^ 2;
t114 = t101 * t102;
t93 = -t101 * t110 + t103 * t106;
t84 = -t98 * t114 + t93 * t97;
t118 = t70 * t84;
t85 = t97 * t114 + t93 * t98;
t109 = t104 * t106;
t92 = t101 * t109 + t103 * t105;
t99 = pkin(12) + qJ(6);
t95 = sin(t99);
t96 = cos(t99);
t76 = t85 * t96 + t92 * t95;
t74 = 0.1e1 / t76 ^ 2;
t75 = t85 * t95 - t92 * t96;
t117 = t74 * t75;
t87 = 0.1e1 / t88 ^ 2;
t116 = t81 * t87;
t115 = t92 * t98;
t111 = t102 * t106;
t108 = t75 ^ 2 * t74 + 0.1e1;
t107 = -t77 * t88 - t78 * t81;
t90 = -t101 * t105 + t103 * t109;
t89 = t104 * t97 + t98 * t112;
t86 = 0.1e1 / t88;
t83 = -t97 * t113 + t91 * t98;
t79 = 0.1e1 / (t81 ^ 2 * t87 + 0.1e1);
t73 = 0.1e1 / t76;
t72 = 0.1e1 / t108;
t69 = 0.1e1 / t71;
t68 = 0.1e1 / (t84 ^ 2 * t70 + 0.1e1);
t67 = (t111 * t116 - t86 * t90) * t97 * t79;
t66 = (t89 * t116 - t83 * t86) * t79;
t65 = (t96 * t117 - t73 * t95) * t84 * t72;
t64 = (t85 * t69 - (t107 * t66 - t77 * t83 + t78 * t89) * t118) * t68;
t1 = [0, t67, t66, t66, 0, 0; 0 (-t92 * t97 * t69 - ((t78 * t111 - t77 * t90) * t97 + t107 * t67) * t118) * t68, t64, t64, 0, 0; 0 ((-t95 * t115 - t93 * t96) * t73 - (-t96 * t115 + t93 * t95) * t117) * t72, t65, t65, 0, t108 * t72;];
Ja_rot  = t1;
