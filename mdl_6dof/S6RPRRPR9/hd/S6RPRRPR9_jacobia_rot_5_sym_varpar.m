% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR9_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:28
% EndTime: 2019-02-26 21:05:29
% DurationCPUTime: 0.24s
% Computational Cost: add. (607->36), mult. (1669->81), div. (55->9), fcn. (2293->15), ass. (0->56)
t101 = cos(pkin(7));
t105 = cos(qJ(1));
t99 = sin(pkin(6));
t117 = t105 * t99;
t102 = cos(pkin(6));
t124 = sin(qJ(1));
t97 = sin(pkin(12));
t111 = t124 * t97;
t100 = cos(pkin(12));
t115 = t105 * t100;
t89 = -t102 * t115 + t111;
t98 = sin(pkin(7));
t125 = t101 * t89 + t98 * t117;
t103 = sin(qJ(3));
t104 = cos(qJ(3));
t110 = t124 * t100;
t118 = t105 * t97;
t90 = t102 * t118 + t110;
t75 = t90 * t103 + t104 * t125;
t107 = t102 * t110 + t118;
t112 = t99 * t124;
t126 = t107 * t101 - t98 * t112;
t91 = -t102 * t111 + t115;
t80 = -t103 * t126 + t91 * t104;
t86 = t101 * t112 + t107 * t98;
t96 = qJ(4) + pkin(13);
t94 = sin(t96);
t95 = cos(t96);
t70 = t80 * t95 + t86 * t94;
t68 = 0.1e1 / t70 ^ 2;
t69 = t80 * t94 - t86 * t95;
t123 = t68 * t69;
t108 = t100 * t101 * t99 + t102 * t98;
t120 = t97 * t99;
t83 = t103 * t120 - t108 * t104;
t74 = atan2(-t75, t83);
t72 = cos(t74);
t122 = t72 * t75;
t71 = sin(t74);
t65 = -t71 * t75 + t72 * t83;
t64 = 0.1e1 / t65 ^ 2;
t79 = t91 * t103 + t104 * t126;
t121 = t79 ^ 2 * t64;
t113 = t69 ^ 2 * t68 + 0.1e1;
t78 = t125 * t103 - t90 * t104;
t85 = t101 * t117 - t89 * t98;
t84 = t108 * t103 + t104 * t120;
t82 = 0.1e1 / t83 ^ 2;
t81 = 0.1e1 / t83;
t73 = 0.1e1 / (t75 ^ 2 * t82 + 0.1e1);
t67 = 0.1e1 / t70;
t66 = 0.1e1 / t113;
t63 = 0.1e1 / t65;
t62 = 0.1e1 / (0.1e1 + t121);
t61 = (t75 * t82 * t84 + t78 * t81) * t73;
t1 = [-t79 * t81 * t73, 0, t61, 0, 0, 0; (-t75 * t63 - (-t71 + (t81 * t122 + t71) * t73) * t121) * t62, 0 (t80 * t63 - (t71 * t78 + t72 * t84 + (-t71 * t83 - t122) * t61) * t79 * t64) * t62, 0, 0, 0; ((t78 * t94 - t85 * t95) * t67 - (t78 * t95 + t85 * t94) * t123) * t66, 0 (t95 * t123 - t94 * t67) * t79 * t66, t113 * t66, 0, 0;];
Ja_rot  = t1;
