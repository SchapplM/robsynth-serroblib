% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRRRR10
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR10_jacobia_rot_3_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_rot_3_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_rot_3_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:17
% EndTime: 2018-11-23 11:27:17
% DurationCPUTime: 0.19s
% Computational Cost: add. (972->44), mult. (1129->86), div. (56->9), fcn. (1288->21), ass. (0->59)
t98 = pkin(7) + qJ(3);
t124 = sin(t98) / 0.2e1;
t105 = cos(qJ(3));
t101 = sin(pkin(6));
t104 = sin(qJ(1));
t118 = t101 * t104;
t103 = sin(qJ(2));
t107 = cos(qJ(1));
t99 = pkin(6) + qJ(2);
t95 = cos(t99) / 0.2e1;
t116 = pkin(6) - qJ(2);
t97 = cos(t116);
t88 = t97 / 0.2e1 + t95;
t79 = -t107 * t103 - t104 * t88;
t106 = cos(qJ(2));
t110 = sin(t116);
t113 = sin(t99) / 0.2e1;
t84 = t113 - t110 / 0.2e1;
t81 = t104 * t84 - t107 * t106;
t115 = pkin(7) - qJ(3);
t109 = sin(t115);
t83 = t124 - t109 / 0.2e1;
t111 = cos(t115);
t114 = cos(t98) / 0.2e1;
t85 = t114 - t111 / 0.2e1;
t64 = -t105 * t81 - t85 * t118 + t79 * t83;
t62 = 0.1e1 / t64 ^ 2;
t102 = sin(qJ(3));
t82 = t124 + t109 / 0.2e1;
t86 = t111 / 0.2e1 + t114;
t63 = -t102 * t81 - t82 * t118 - t79 * t86;
t123 = t62 * t63;
t122 = t63 ^ 2 * t62;
t100 = sin(pkin(7));
t119 = cos(pkin(7));
t112 = t101 * t119;
t76 = t104 * t103 - t107 * t88;
t69 = -t76 * t100 + t107 * t112;
t75 = -(t113 + t110 / 0.2e1) * t100 + cos(pkin(6)) * t119;
t68 = atan2(t69, t75);
t66 = cos(t68);
t121 = t66 * t69;
t65 = sin(t68);
t59 = t65 * t69 + t66 * t75;
t58 = 0.1e1 / t59 ^ 2;
t70 = t79 * t100 - t104 * t112;
t120 = t70 ^ 2 * t58;
t117 = t101 * t107;
t108 = -t104 * t106 - t107 * t84;
t87 = t95 - t97 / 0.2e1;
t74 = 0.1e1 / t75 ^ 2;
t73 = 0.1e1 / t75;
t67 = 0.1e1 / (t69 ^ 2 * t74 + 0.1e1);
t61 = 0.1e1 / t64;
t60 = 0.1e1 / (0.1e1 + t122);
t57 = 0.1e1 / t59;
t56 = 0.1e1 / (0.1e1 + t120);
t55 = (t69 * t74 * t87 + t108 * t73) * t67 * t100;
t1 = [t70 * t73 * t67, t55, 0, 0, 0, 0; (t69 * t57 + (t65 + (t73 * t121 - t65) * t67) * t120) * t56 (-t81 * t100 * t57 + ((-t65 * t75 + t121) * t55 + (t108 * t65 - t66 * t87) * t100) * t70 * t58) * t56, 0, 0, 0, 0; ((t102 * t108 - t82 * t117 - t76 * t86) * t61 - (t105 * t108 - t85 * t117 + t76 * t83) * t123) * t60 ((t79 * t102 - t81 * t86) * t61 - (t79 * t105 + t81 * t83) * t123) * t60 (t61 * t64 + t122) * t60, 0, 0, 0;];
Ja_rot  = t1;
