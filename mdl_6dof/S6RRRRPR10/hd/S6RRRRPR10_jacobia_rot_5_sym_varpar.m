% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR10_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:48
% EndTime: 2019-02-26 22:35:48
% DurationCPUTime: 0.21s
% Computational Cost: add. (992->33), mult. (1383->75), div. (104->9), fcn. (2002->11), ass. (0->51)
t88 = cos(pkin(6));
t89 = sin(qJ(2));
t92 = cos(qJ(1));
t95 = t92 * t89;
t90 = sin(qJ(1));
t91 = cos(qJ(2));
t96 = t90 * t91;
t79 = t88 * t95 + t96;
t86 = qJ(3) + qJ(4);
t84 = sin(t86);
t85 = cos(t86);
t87 = sin(pkin(6));
t98 = t87 * t92;
t67 = t79 * t84 + t85 * t98;
t101 = t87 * t89;
t74 = t84 * t101 - t88 * t85;
t65 = atan2(-t67, t74);
t62 = sin(t65);
t63 = cos(t65);
t60 = -t62 * t67 + t63 * t74;
t59 = 0.1e1 / t60 ^ 2;
t100 = t87 * t90;
t94 = t92 * t91;
t97 = t90 * t89;
t81 = -t88 * t97 + t94;
t70 = -t85 * t100 + t81 * t84;
t106 = t59 * t70;
t105 = t63 * t67;
t73 = 0.1e1 / t74 ^ 2;
t104 = t67 * t73;
t103 = t70 ^ 2 * t59;
t71 = t84 * t100 + t81 * t85;
t80 = t88 * t96 + t95;
t77 = 0.1e1 / t80 ^ 2;
t102 = t71 * t77;
t99 = t87 * t91;
t69 = t79 * t85 - t84 * t98;
t93 = -t62 * t74 - t105;
t78 = t88 * t94 - t97;
t76 = 0.1e1 / t80;
t75 = t85 * t101 + t88 * t84;
t72 = 0.1e1 / t74;
t66 = 0.1e1 / (t71 ^ 2 * t77 + 0.1e1);
t64 = 0.1e1 / (t67 ^ 2 * t73 + 0.1e1);
t61 = t70 * t76 * t66;
t58 = 0.1e1 / t60;
t57 = 0.1e1 / (0.1e1 + t103);
t56 = (t99 * t104 - t72 * t78) * t84 * t64;
t55 = (t75 * t104 - t69 * t72) * t64;
t54 = (t71 * t58 - (t93 * t55 - t62 * t69 + t63 * t75) * t106) * t57;
t1 = [-t70 * t72 * t64, t56, t55, t55, 0, 0; (-t67 * t58 - (-t62 + (t72 * t105 + t62) * t64) * t103) * t57 (-t80 * t84 * t58 - ((-t62 * t78 + t63 * t99) * t84 + t93 * t56) * t106) * t57, t54, t54, 0, 0; (-t78 * t102 - t69 * t76) * t66 (-t76 * t80 * t85 - t81 * t102) * t66, -t61, -t61, 0, 0;];
Ja_rot  = t1;
