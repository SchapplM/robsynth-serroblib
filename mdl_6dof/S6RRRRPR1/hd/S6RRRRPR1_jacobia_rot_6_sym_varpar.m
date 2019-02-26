% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR1
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
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:46
% EndTime: 2019-02-26 22:30:46
% DurationCPUTime: 0.09s
% Computational Cost: add. (1004->21), mult. (440->54), div. (106->9), fcn. (646->9), ass. (0->38)
t80 = qJ(2) + qJ(3) + qJ(4) + pkin(11);
t79 = cos(t80);
t78 = sin(t80);
t82 = sin(qJ(1));
t90 = t82 * t78;
t69 = atan2(-t90, -t79);
t67 = sin(t69);
t68 = cos(t69);
t64 = -t67 * t90 - t68 * t79;
t63 = 0.1e1 / t64 ^ 2;
t84 = cos(qJ(1));
t96 = t63 * t84 ^ 2;
t95 = t67 * t79;
t83 = cos(qJ(6));
t86 = t84 * t83;
t81 = sin(qJ(6));
t89 = t82 * t81;
t74 = t79 * t86 + t89;
t71 = 0.1e1 / t74 ^ 2;
t87 = t84 * t81;
t88 = t82 * t83;
t73 = t79 * t87 - t88;
t94 = t71 * t73;
t75 = t78 ^ 2;
t93 = t75 / t79 ^ 2;
t92 = t78 * t84;
t72 = 0.1e1 / (t82 ^ 2 * t93 + 0.1e1);
t91 = t82 * t72;
t85 = t73 ^ 2 * t71 + 0.1e1;
t76 = 0.1e1 / t79;
t70 = 0.1e1 / t74;
t66 = 0.1e1 / t85;
t65 = (0.1e1 + t93) * t91;
t62 = 0.1e1 / t64;
t61 = 0.1e1 / (t75 * t96 + 0.1e1);
t60 = (-t70 * t81 + t83 * t94) * t66 * t92;
t59 = (t79 * t62 - (-t82 * t95 + t68 * t78 + (-t68 * t90 + t95) * t65) * t78 * t63) * t84 * t61;
t1 = [t76 * t72 * t92, t65, t65, t65, 0, 0; (-t62 * t90 - (-t68 * t75 * t76 * t91 + (t72 - 0.1e1) * t78 * t67) * t78 * t96) * t61, t59, t59, t59, 0, 0; ((-t79 * t89 - t86) * t70 - (-t79 * t88 + t87) * t94) * t66, t60, t60, t60, 0, t85 * t66;];
Ja_rot  = t1;
