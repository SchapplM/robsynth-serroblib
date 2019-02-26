% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:32:17
% EndTime: 2019-02-26 22:32:17
% DurationCPUTime: 0.08s
% Computational Cost: add. (531->22), mult. (359->54), div. (84->9), fcn. (524->9), ass. (0->40)
t82 = qJ(2) + qJ(3);
t81 = cos(t82);
t80 = sin(t82);
t83 = sin(qJ(1));
t88 = t83 * t80;
t72 = atan2(-t88, -t81);
t70 = sin(t72);
t71 = cos(t72);
t64 = -t70 * t88 - t71 * t81;
t63 = 0.1e1 / t64 ^ 2;
t84 = cos(qJ(1));
t96 = t63 * t84 ^ 2;
t79 = qJ(4) + pkin(11) + qJ(6);
t75 = cos(t79);
t86 = t84 * t75;
t74 = sin(t79);
t90 = t83 * t74;
t69 = t81 * t86 + t90;
t67 = 0.1e1 / t69 ^ 2;
t87 = t84 * t74;
t89 = t83 * t75;
t68 = t81 * t87 - t89;
t95 = t67 * t68;
t94 = t70 * t81;
t76 = t80 ^ 2;
t93 = t76 / t81 ^ 2;
t92 = t80 * t84;
t73 = 0.1e1 / (t83 ^ 2 * t93 + 0.1e1);
t91 = t83 * t73;
t85 = t68 ^ 2 * t67 + 0.1e1;
t77 = 0.1e1 / t81;
t66 = 0.1e1 / t69;
t65 = (0.1e1 + t93) * t91;
t62 = 0.1e1 / t64;
t61 = 0.1e1 / t85;
t60 = 0.1e1 / (t76 * t96 + 0.1e1);
t59 = t85 * t61;
t58 = (-t66 * t74 + t75 * t95) * t61 * t92;
t57 = (t81 * t62 - (-t83 * t94 + t71 * t80 + (-t71 * t88 + t94) * t65) * t80 * t63) * t84 * t60;
t1 = [t77 * t73 * t92, t65, t65, 0, 0, 0; (-t62 * t88 - (-t71 * t76 * t77 * t91 + (t73 - 0.1e1) * t80 * t70) * t80 * t96) * t60, t57, t57, 0, 0, 0; ((-t81 * t90 - t86) * t66 - (-t81 * t89 + t87) * t95) * t61, t58, t58, t59, 0, t59;];
Ja_rot  = t1;
