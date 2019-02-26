% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP3
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

function Ja_rot = S6RRRRRP3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:40:41
% EndTime: 2019-02-26 22:40:42
% DurationCPUTime: 0.09s
% Computational Cost: add. (453->22), mult. (359->54), div. (84->9), fcn. (524->9), ass. (0->40)
t81 = qJ(2) + qJ(3);
t79 = cos(t81);
t77 = sin(t81);
t82 = sin(qJ(1));
t88 = t82 * t77;
t71 = atan2(-t88, -t79);
t69 = sin(t71);
t70 = cos(t71);
t62 = -t69 * t88 - t70 * t79;
t61 = 0.1e1 / t62 ^ 2;
t83 = cos(qJ(1));
t95 = t61 * t83 ^ 2;
t80 = qJ(4) + qJ(5);
t78 = cos(t80);
t85 = t83 * t78;
t76 = sin(t80);
t89 = t82 * t76;
t68 = t79 * t85 + t89;
t66 = 0.1e1 / t68 ^ 2;
t86 = t83 * t76;
t87 = t82 * t78;
t67 = t79 * t86 - t87;
t94 = t66 * t67;
t93 = t69 * t79;
t73 = t77 ^ 2;
t92 = t73 / t79 ^ 2;
t91 = t77 * t83;
t72 = 0.1e1 / (t82 ^ 2 * t92 + 0.1e1);
t90 = t82 * t72;
t84 = t67 ^ 2 * t66 + 0.1e1;
t74 = 0.1e1 / t79;
t65 = 0.1e1 / t68;
t64 = (0.1e1 + t92) * t90;
t63 = 0.1e1 / t84;
t60 = 0.1e1 / t62;
t59 = 0.1e1 / (t73 * t95 + 0.1e1);
t58 = t84 * t63;
t57 = (-t65 * t76 + t78 * t94) * t63 * t91;
t56 = (t79 * t60 - (-t82 * t93 + t70 * t77 + (-t70 * t88 + t93) * t64) * t77 * t61) * t83 * t59;
t1 = [t74 * t72 * t91, t64, t64, 0, 0, 0; (-t60 * t88 - (-t70 * t73 * t74 * t90 + (t72 - 0.1e1) * t77 * t69) * t77 * t95) * t59, t56, t56, 0, 0, 0; ((-t79 * t89 - t85) * t65 - (-t79 * t87 + t86) * t94) * t63, t57, t57, t58, t58, 0;];
Ja_rot  = t1;
