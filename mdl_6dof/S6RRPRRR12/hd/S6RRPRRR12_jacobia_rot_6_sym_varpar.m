% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR12_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:00:24
% EndTime: 2019-02-26 22:00:24
% DurationCPUTime: 0.22s
% Computational Cost: add. (1185->37), mult. (1738->90), div. (115->9), fcn. (2485->13), ass. (0->58)
t103 = cos(qJ(1));
t96 = sin(pkin(6));
t112 = t103 * t96;
t102 = cos(qJ(2));
t108 = t103 * t102;
t100 = sin(qJ(1));
t99 = sin(qJ(2));
t114 = t100 * t99;
t97 = cos(pkin(6));
t89 = -t97 * t108 + t114;
t95 = qJ(4) + qJ(5);
t93 = sin(t95);
t94 = cos(t95);
t81 = t93 * t112 + t89 * t94;
t113 = t102 * t96;
t87 = t94 * t113 + t97 * t93;
t78 = atan2(t81, t87);
t73 = sin(t78);
t74 = cos(t78);
t69 = t73 * t81 + t74 * t87;
t68 = 0.1e1 / t69 ^ 2;
t109 = t100 * t102;
t111 = t103 * t99;
t104 = t97 * t109 + t111;
t115 = t100 * t96;
t79 = -t104 * t94 + t93 * t115;
t122 = t68 * t79;
t101 = cos(qJ(6));
t91 = -t97 * t114 + t108;
t98 = sin(qJ(6));
t117 = t91 * t98;
t80 = t104 * t93 + t94 * t115;
t76 = t80 * t101 + t117;
t72 = 0.1e1 / t76 ^ 2;
t110 = t91 * t101;
t75 = t80 * t98 - t110;
t121 = t72 * t75;
t120 = t74 * t81;
t119 = t79 ^ 2 * t68;
t86 = 0.1e1 / t87 ^ 2;
t118 = t81 * t86;
t116 = t96 * t99;
t107 = t75 ^ 2 * t72 + 0.1e1;
t106 = -t73 * t87 + t120;
t105 = t94 * t112 - t89 * t93;
t90 = t97 * t111 + t109;
t88 = -t93 * t113 + t97 * t94;
t85 = 0.1e1 / t87;
t77 = 0.1e1 / (t81 ^ 2 * t86 + 0.1e1);
t71 = 0.1e1 / t76;
t70 = 0.1e1 / t107;
t67 = 0.1e1 / t69;
t66 = 0.1e1 / (0.1e1 + t119);
t65 = (t116 * t118 + t85 * t90) * t94 * t77;
t64 = (t105 * t85 - t88 * t118) * t77;
t63 = (t101 * t121 - t98 * t71) * t79 * t70;
t62 = (t80 * t67 - (t105 * t73 + t106 * t64 + t74 * t88) * t122) * t66;
t1 = [-t79 * t85 * t77, t65, 0, t64, t64, 0; (t81 * t67 - (-t73 + (-t85 * t120 + t73) * t77) * t119) * t66 (-t91 * t94 * t67 - ((-t74 * t116 + t73 * t90) * t94 + t106 * t65) * t122) * t66, 0, t62, t62, 0; ((t90 * t101 + t105 * t98) * t71 - (t101 * t105 - t90 * t98) * t121) * t70 ((t104 * t101 + t93 * t117) * t71 - (-t104 * t98 + t93 * t110) * t121) * t70, 0, t63, t63, t107 * t70;];
Ja_rot  = t1;
