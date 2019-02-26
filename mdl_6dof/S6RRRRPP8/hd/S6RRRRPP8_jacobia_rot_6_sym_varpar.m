% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:29:34
% EndTime: 2019-02-26 22:29:34
% DurationCPUTime: 0.18s
% Computational Cost: add. (463->36), mult. (1311->91), div. (86->9), fcn. (1874->13), ass. (0->54)
t119 = sin(qJ(1));
t100 = cos(qJ(3));
t102 = cos(qJ(1));
t94 = sin(pkin(6));
t109 = t102 * t94;
t101 = cos(qJ(2));
t104 = t119 * t101;
t98 = sin(qJ(2));
t108 = t102 * t98;
t95 = cos(pkin(6));
t89 = t95 * t108 + t104;
t97 = sin(qJ(3));
t78 = t100 * t109 + t89 * t97;
t112 = t94 * t98;
t86 = t95 * t100 - t97 * t112;
t77 = atan2(t78, t86);
t74 = sin(t77);
t75 = cos(t77);
t67 = t74 * t78 + t75 * t86;
t66 = 0.1e1 / t67 ^ 2;
t106 = t94 * t119;
t105 = t119 * t98;
t107 = t102 * t101;
t91 = -t95 * t105 + t107;
t81 = -t100 * t106 + t91 * t97;
t118 = t66 * t81;
t83 = t91 * t100 + t97 * t106;
t90 = t95 * t104 + t108;
t96 = sin(qJ(4));
t99 = cos(qJ(4));
t73 = t83 * t99 + t90 * t96;
t71 = 0.1e1 / t73 ^ 2;
t72 = -t83 * t96 + t90 * t99;
t117 = t72 ^ 2 * t71;
t116 = t71 * t72;
t115 = t75 * t78;
t85 = 0.1e1 / t86 ^ 2;
t114 = t78 * t85;
t113 = t81 ^ 2 * t66;
t111 = t100 * t90;
t110 = t101 * t94;
t79 = t89 * t100 - t97 * t109;
t103 = -t74 * t86 + t115;
t88 = t95 * t107 - t105;
t87 = -t100 * t112 - t95 * t97;
t84 = 0.1e1 / t86;
t76 = 0.1e1 / (t78 ^ 2 * t85 + 0.1e1);
t70 = 0.1e1 / t73;
t68 = 0.1e1 / (0.1e1 + t117);
t65 = 0.1e1 / t67;
t64 = 0.1e1 / (0.1e1 + t113);
t63 = (t110 * t114 + t84 * t88) * t97 * t76;
t62 = (-t87 * t114 + t79 * t84) * t76;
t1 = [t81 * t84 * t76, t63, t62, 0, 0, 0; (t78 * t65 + (t74 + (t84 * t115 - t74) * t76) * t113) * t64 (t90 * t97 * t65 + ((-t75 * t110 + t74 * t88) * t97 + t103 * t63) * t118) * t64 (-t83 * t65 + (t103 * t62 + t74 * t79 + t75 * t87) * t118) * t64, 0, 0, 0; ((t79 * t96 + t88 * t99) * t70 - (-t79 * t99 + t88 * t96) * t116) * t68 ((t96 * t111 + t91 * t99) * t70 - (-t99 * t111 + t91 * t96) * t116) * t68 (t99 * t116 + t96 * t70) * t81 * t68 (-t70 * t73 - t117) * t68, 0, 0;];
Ja_rot  = t1;
