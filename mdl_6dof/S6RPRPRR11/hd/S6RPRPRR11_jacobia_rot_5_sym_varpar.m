% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR11_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:31
% EndTime: 2019-02-26 20:54:31
% DurationCPUTime: 0.24s
% Computational Cost: add. (607->36), mult. (1669->83), div. (55->9), fcn. (2293->15), ass. (0->54)
t100 = cos(pkin(7));
t104 = cos(qJ(1));
t98 = sin(pkin(6));
t115 = t104 * t98;
t101 = cos(pkin(6));
t113 = t101 * t104;
t121 = sin(qJ(1));
t96 = sin(pkin(12));
t99 = cos(pkin(12));
t88 = -t99 * t113 + t121 * t96;
t97 = sin(pkin(7));
t122 = t100 * t88 + t97 * t115;
t102 = sin(qJ(3));
t103 = cos(qJ(3));
t89 = t96 * t113 + t121 * t99;
t74 = t89 * t102 + t103 * t122;
t109 = t101 * t121;
t106 = t104 * t96 + t99 * t109;
t110 = t98 * t121;
t123 = t106 * t100 - t97 * t110;
t90 = t104 * t99 - t96 * t109;
t79 = -t102 * t123 + t90 * t103;
t85 = t100 * t110 + t106 * t97;
t95 = pkin(13) + qJ(5);
t93 = sin(t95);
t94 = cos(t95);
t69 = t79 * t94 + t85 * t93;
t67 = 0.1e1 / t69 ^ 2;
t68 = t79 * t93 - t85 * t94;
t120 = t67 * t68;
t107 = t100 * t98 * t99 + t101 * t97;
t117 = t96 * t98;
t82 = t102 * t117 - t107 * t103;
t73 = atan2(-t74, t82);
t71 = cos(t73);
t119 = t71 * t74;
t70 = sin(t73);
t64 = -t70 * t74 + t71 * t82;
t63 = 0.1e1 / t64 ^ 2;
t78 = t90 * t102 + t103 * t123;
t118 = t78 ^ 2 * t63;
t111 = t68 ^ 2 * t67 + 0.1e1;
t77 = t122 * t102 - t89 * t103;
t84 = t100 * t115 - t88 * t97;
t83 = t107 * t102 + t103 * t117;
t81 = 0.1e1 / t82 ^ 2;
t80 = 0.1e1 / t82;
t72 = 0.1e1 / (t74 ^ 2 * t81 + 0.1e1);
t66 = 0.1e1 / t69;
t65 = 0.1e1 / t111;
t62 = 0.1e1 / t64;
t61 = 0.1e1 / (0.1e1 + t118);
t60 = (t74 * t81 * t83 + t77 * t80) * t72;
t1 = [-t78 * t80 * t72, 0, t60, 0, 0, 0; (-t74 * t62 - (-t70 + (t80 * t119 + t70) * t72) * t118) * t61, 0 (t79 * t62 - (t70 * t77 + t71 * t83 + (-t70 * t82 - t119) * t60) * t78 * t63) * t61, 0, 0, 0; ((t77 * t93 - t84 * t94) * t66 - (t77 * t94 + t84 * t93) * t120) * t65, 0 (t94 * t120 - t93 * t66) * t78 * t65, 0, t111 * t65, 0;];
Ja_rot  = t1;
