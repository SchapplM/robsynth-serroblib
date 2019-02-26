% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP2
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
% Datum: 2019-02-26 22:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRP2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:40:10
% EndTime: 2019-02-26 22:40:10
% DurationCPUTime: 0.13s
% Computational Cost: add. (933->27), mult. (867->69), div. (172->11), fcn. (1298->9), ass. (0->43)
t85 = qJ(2) + qJ(3) + qJ(4);
t83 = sin(t85);
t105 = t83 ^ 2;
t84 = cos(t85);
t91 = cos(qJ(5));
t92 = cos(qJ(1));
t94 = t92 * t91;
t89 = sin(qJ(5));
t90 = sin(qJ(1));
t97 = t90 * t89;
t73 = t84 * t97 + t94;
t99 = t83 * t89;
t70 = atan2(-t73, t99);
t66 = sin(t70);
t67 = cos(t70);
t65 = -t66 * t73 + t67 * t99;
t64 = 0.1e1 / t65 ^ 2;
t95 = t92 * t89;
t96 = t90 * t91;
t76 = t84 * t95 - t96;
t104 = t64 * t76;
t102 = t67 * t73;
t101 = t76 ^ 2 * t64;
t79 = 0.1e1 / t83;
t86 = 0.1e1 / t89;
t100 = t79 * t86;
t98 = t83 * t92;
t77 = t84 * t94 + t97;
t72 = 0.1e1 / t77 ^ 2;
t93 = t92 ^ 2 * t105 * t72;
t87 = 0.1e1 / t89 ^ 2;
t80 = 0.1e1 / t105;
t75 = t84 * t96 - t95;
t71 = 0.1e1 / t77;
t69 = 0.1e1 / (t73 ^ 2 * t80 * t87 + 0.1e1);
t68 = 0.1e1 / (0.1e1 + t93);
t63 = 0.1e1 / t65;
t62 = (t73 * t80 * t84 * t86 + t90) * t69;
t61 = 0.1e1 / (0.1e1 + t101);
t60 = (t73 * t87 * t91 - t75 * t86) * t79 * t69;
t59 = (-t71 * t84 * t92 - t91 * t93) * t68;
t58 = (t62 * t102 * t104 + (-t63 * t98 - (t67 * t84 + (-t62 + t90) * t83 * t66) * t104) * t89) * t61;
t1 = [-t76 * t69 * t100, t62, t62, t62, t60, 0; (-t73 * t63 - (-t66 + (t100 * t102 + t66) * t69) * t101) * t61, t58, t58, t58 (t77 * t63 - (t67 * t83 * t91 - t66 * t75 + (-t66 * t99 - t102) * t60) * t104) * t61, 0; (-t72 * t75 * t92 + t71 * t90) * t83 * t68, t59, t59, t59, -t76 * t72 * t68 * t98, 0;];
Ja_rot  = t1;
