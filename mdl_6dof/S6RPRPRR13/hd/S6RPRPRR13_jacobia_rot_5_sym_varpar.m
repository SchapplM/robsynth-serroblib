% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR13_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:38
% EndTime: 2019-02-26 20:55:38
% DurationCPUTime: 0.24s
% Computational Cost: add. (559->36), mult. (1669->85), div. (55->9), fcn. (2293->15), ass. (0->55)
t90 = sin(pkin(6));
t99 = cos(qJ(1));
t110 = t90 * t99;
t89 = sin(pkin(7));
t104 = t89 * t110;
t91 = cos(pkin(12));
t105 = t99 * t91;
t88 = sin(pkin(12));
t96 = sin(qJ(1));
t108 = t96 * t88;
t93 = cos(pkin(6));
t82 = -t93 * t105 + t108;
t106 = t99 * t88;
t107 = t96 * t91;
t83 = t93 * t106 + t107;
t92 = cos(pkin(7));
t95 = sin(qJ(3));
t98 = cos(qJ(3));
t70 = (t82 * t92 + t104) * t98 + t83 * t95;
t111 = t90 * t96;
t84 = -t93 * t107 - t106;
t101 = t89 * t111 + t84 * t92;
t85 = -t93 * t108 + t105;
t74 = -t101 * t98 + t85 * t95;
t81 = t92 * t111 - t84 * t89;
t94 = sin(qJ(5));
t97 = cos(qJ(5));
t65 = t74 * t94 + t81 * t97;
t63 = 0.1e1 / t65 ^ 2;
t64 = -t74 * t97 + t81 * t94;
t116 = t63 * t64;
t109 = t92 * t95;
t100 = t95 * t104 + t82 * t109 - t83 * t98;
t112 = t89 * t93;
t79 = t95 * t112 + (t91 * t109 + t88 * t98) * t90;
t69 = atan2(t100, t79);
t67 = cos(t69);
t115 = t67 * t100;
t66 = sin(t69);
t60 = t100 * t66 + t67 * t79;
t59 = 0.1e1 / t60 ^ 2;
t75 = t101 * t95 + t85 * t98;
t114 = t75 ^ 2 * t59;
t103 = t64 ^ 2 * t63 + 0.1e1;
t80 = t92 * t110 - t82 * t89;
t78 = t98 * t112 + (t91 * t92 * t98 - t88 * t95) * t90;
t77 = 0.1e1 / t79 ^ 2;
t76 = 0.1e1 / t79;
t68 = 0.1e1 / (t100 ^ 2 * t77 + 0.1e1);
t62 = 0.1e1 / t65;
t61 = 0.1e1 / t103;
t58 = 0.1e1 / t60;
t57 = 0.1e1 / (0.1e1 + t114);
t56 = (-t100 * t77 * t78 + t70 * t76) * t68;
t1 = [-t75 * t76 * t68, 0, t56, 0, 0, 0; (t100 * t58 - (-t66 + (-t76 * t115 + t66) * t68) * t114) * t57, 0 (-t74 * t58 - (t66 * t70 + t67 * t78 + (-t66 * t79 + t115) * t56) * t75 * t59) * t57, 0, 0, 0; ((t70 * t97 + t80 * t94) * t62 - (-t70 * t94 + t80 * t97) * t116) * t61, 0 (-t94 * t116 - t97 * t62) * t75 * t61, 0, t103 * t61, 0;];
Ja_rot  = t1;
