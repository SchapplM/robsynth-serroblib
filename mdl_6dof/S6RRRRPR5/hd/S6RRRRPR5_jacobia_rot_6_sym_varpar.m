% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:33:03
% EndTime: 2019-02-26 22:33:03
% DurationCPUTime: 0.16s
% Computational Cost: add. (434->26), mult. (576->67), div. (85->9), fcn. (827->11), ass. (0->45)
t100 = cos(qJ(1));
t96 = sin(qJ(4));
t104 = t100 * t96;
t97 = sin(qJ(1));
t99 = cos(qJ(4));
t106 = t97 * t99;
t94 = qJ(2) + qJ(3);
t93 = cos(t94);
t82 = t93 * t104 - t106;
t103 = t100 * t99;
t107 = t97 * t96;
t83 = t93 * t103 + t107;
t95 = sin(qJ(6));
t98 = cos(qJ(6));
t79 = t82 * t95 + t83 * t98;
t77 = 0.1e1 / t79 ^ 2;
t78 = -t82 * t98 + t83 * t95;
t113 = t77 * t78;
t112 = t78 ^ 2 * t77;
t92 = sin(t94);
t108 = t97 * t92;
t87 = atan2(t108, t93);
t84 = sin(t87);
t111 = t84 * t93;
t89 = t92 ^ 2;
t110 = t89 / t93 ^ 2;
t86 = 0.1e1 / (t97 ^ 2 * t110 + 0.1e1);
t109 = t97 * t86;
t105 = t100 * t92;
t85 = cos(t87);
t74 = t84 * t108 + t85 * t93;
t73 = 0.1e1 / t74 ^ 2;
t102 = t73 * t100 ^ 2;
t101 = 0.1e1 + t112;
t90 = 0.1e1 / t93;
t81 = -t93 * t106 + t104;
t80 = -t93 * t107 - t103;
t76 = 0.1e1 / t79;
t75 = (0.1e1 + t110) * t109;
t72 = 0.1e1 / t74;
t71 = 0.1e1 / (t89 * t102 + 0.1e1);
t70 = 0.1e1 / t101;
t69 = ((-t95 * t99 + t96 * t98) * t76 - (-t95 * t96 - t98 * t99) * t113) * t70 * t105;
t68 = (-t93 * t72 + (t97 * t111 - t85 * t92 + (t85 * t108 - t111) * t75) * t92 * t73) * t71 * t100;
t1 = [t90 * t86 * t105, t75, t75, 0, 0, 0; (t72 * t108 + (t85 * t89 * t90 * t109 + (-t86 + 0.1e1) * t92 * t84) * t92 * t102) * t71, t68, t68, 0, 0, 0; ((-t80 * t98 + t81 * t95) * t76 - (t80 * t95 + t81 * t98) * t113) * t70, t69, t69 (-t79 * t76 - t112) * t70, 0, t101 * t70;];
Ja_rot  = t1;
