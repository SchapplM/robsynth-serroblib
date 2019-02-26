% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR12_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:22:35
% EndTime: 2019-02-26 22:22:35
% DurationCPUTime: 0.20s
% Computational Cost: add. (650->38), mult. (1383->91), div. (90->9), fcn. (1970->13), ass. (0->56)
t92 = sin(pkin(6));
t99 = cos(qJ(1));
t106 = t92 * t99;
t95 = sin(qJ(2));
t103 = t99 * t95;
t96 = sin(qJ(1));
t98 = cos(qJ(2));
t104 = t96 * t98;
t93 = cos(pkin(6));
t84 = t93 * t103 + t104;
t94 = sin(qJ(3));
t97 = cos(qJ(3));
t73 = t97 * t106 + t84 * t94;
t109 = t92 * t94;
t81 = t95 * t109 - t93 * t97;
t72 = atan2(-t73, t81);
t69 = sin(t72);
t70 = cos(t72);
t63 = -t69 * t73 + t70 * t81;
t62 = 0.1e1 / t63 ^ 2;
t108 = t92 * t97;
t102 = t99 * t98;
t105 = t96 * t95;
t86 = -t93 * t105 + t102;
t77 = -t96 * t108 + t86 * t94;
t115 = t62 * t77;
t78 = t96 * t109 + t86 * t97;
t85 = t93 * t104 + t103;
t91 = pkin(12) + qJ(5) + qJ(6);
t89 = sin(t91);
t90 = cos(t91);
t68 = t78 * t90 + t85 * t89;
t66 = 0.1e1 / t68 ^ 2;
t67 = t78 * t89 - t85 * t90;
t114 = t66 * t67;
t113 = t70 * t73;
t80 = 0.1e1 / t81 ^ 2;
t112 = t73 * t80;
t111 = t77 ^ 2 * t62;
t110 = t85 * t97;
t107 = t92 * t98;
t101 = t67 ^ 2 * t66 + 0.1e1;
t75 = -t94 * t106 + t84 * t97;
t100 = -t69 * t81 - t113;
t83 = t93 * t102 - t105;
t82 = t95 * t108 + t93 * t94;
t79 = 0.1e1 / t81;
t71 = 0.1e1 / (t73 ^ 2 * t80 + 0.1e1);
t65 = 0.1e1 / t68;
t64 = 0.1e1 / t101;
t61 = 0.1e1 / t63;
t60 = 0.1e1 / (0.1e1 + t111);
t59 = (t107 * t112 - t79 * t83) * t94 * t71;
t58 = (t82 * t112 - t75 * t79) * t71;
t57 = t101 * t64;
t1 = [-t77 * t79 * t71, t59, t58, 0, 0, 0; (-t73 * t61 - (-t69 + (t79 * t113 + t69) * t71) * t111) * t60 (-t85 * t94 * t61 - ((t70 * t107 - t69 * t83) * t94 + t100 * t59) * t115) * t60 (t78 * t61 - (t100 * t58 - t69 * t75 + t70 * t82) * t115) * t60, 0, 0, 0; ((-t75 * t89 - t83 * t90) * t65 - (-t75 * t90 + t83 * t89) * t114) * t64 ((-t89 * t110 - t86 * t90) * t65 - (-t90 * t110 + t86 * t89) * t114) * t64 (t90 * t114 - t89 * t65) * t77 * t64, 0, t57, t57;];
Ja_rot  = t1;
