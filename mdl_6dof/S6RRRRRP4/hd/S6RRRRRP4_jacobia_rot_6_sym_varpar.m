% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRP4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:41:10
% EndTime: 2019-02-26 22:41:10
% DurationCPUTime: 0.18s
% Computational Cost: add. (1093->29), mult. (868->67), div. (175->11), fcn. (1310->9), ass. (0->45)
t98 = qJ(2) + qJ(3);
t93 = sin(t98);
t112 = t93 ^ 2;
t97 = qJ(4) + qJ(5);
t92 = sin(t97);
t106 = t93 * t92;
t100 = cos(qJ(1));
t99 = sin(qJ(1));
t105 = t99 * t92;
t94 = cos(t97);
t95 = cos(t98);
t80 = t100 * t94 + t95 * t105;
t76 = atan2(-t80, t106);
t73 = sin(t76);
t74 = cos(t76);
t71 = t74 * t106 - t73 * t80;
t70 = 0.1e1 / t71 ^ 2;
t102 = t100 * t95;
t104 = t99 * t94;
t83 = t92 * t102 - t104;
t111 = t70 * t83;
t109 = t74 * t80;
t108 = t83 ^ 2 * t70;
t87 = 0.1e1 / t92;
t90 = 0.1e1 / t93;
t107 = t87 * t90;
t103 = t100 * t93;
t84 = t94 * t102 + t105;
t79 = 0.1e1 / t84 ^ 2;
t101 = t100 ^ 2 * t112 * t79;
t91 = 0.1e1 / t112;
t88 = 0.1e1 / t92 ^ 2;
t82 = -t100 * t92 + t95 * t104;
t78 = 0.1e1 / t84;
t77 = 0.1e1 / (0.1e1 + t101);
t75 = 0.1e1 / (t80 ^ 2 * t91 * t88 + 0.1e1);
t72 = t83 * t79 * t77 * t103;
t69 = 0.1e1 / t71;
t68 = (t80 * t87 * t91 * t95 + t99) * t75;
t67 = (-t94 * t101 - t78 * t102) * t77;
t66 = 0.1e1 / (0.1e1 + t108);
t65 = (t80 * t88 * t94 - t82 * t87) * t90 * t75;
t64 = (t68 * t109 * t111 + (-t69 * t103 - (t74 * t95 + (-t68 + t99) * t93 * t73) * t111) * t92) * t66;
t63 = (t84 * t69 - (t74 * t93 * t94 - t73 * t82 + (-t73 * t106 - t109) * t65) * t111) * t66;
t1 = [-t83 * t75 * t107, t68, t68, t65, t65, 0; (-t80 * t69 - (-t73 + (t107 * t109 + t73) * t75) * t108) * t66, t64, t64, t63, t63, 0; (-t100 * t79 * t82 + t78 * t99) * t93 * t77, t67, t67, -t72, -t72, 0;];
Ja_rot  = t1;
