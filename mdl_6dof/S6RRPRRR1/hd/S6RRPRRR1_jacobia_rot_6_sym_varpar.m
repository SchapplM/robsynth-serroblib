% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:54:01
% EndTime: 2019-02-26 21:54:01
% DurationCPUTime: 0.09s
% Computational Cost: add. (1004->21), mult. (440->54), div. (106->9), fcn. (646->9), ass. (0->38)
t77 = qJ(2) + pkin(11) + qJ(4) + qJ(5);
t76 = cos(t77);
t75 = sin(t77);
t79 = sin(qJ(1));
t87 = t79 * t75;
t66 = atan2(-t87, -t76);
t64 = sin(t66);
t65 = cos(t66);
t61 = -t64 * t87 - t65 * t76;
t60 = 0.1e1 / t61 ^ 2;
t81 = cos(qJ(1));
t93 = t60 * t81 ^ 2;
t92 = t64 * t76;
t80 = cos(qJ(6));
t83 = t81 * t80;
t78 = sin(qJ(6));
t86 = t79 * t78;
t71 = t76 * t83 + t86;
t68 = 0.1e1 / t71 ^ 2;
t84 = t81 * t78;
t85 = t79 * t80;
t70 = t76 * t84 - t85;
t91 = t68 * t70;
t72 = t75 ^ 2;
t90 = t72 / t76 ^ 2;
t89 = t75 * t81;
t69 = 0.1e1 / (t79 ^ 2 * t90 + 0.1e1);
t88 = t79 * t69;
t82 = t70 ^ 2 * t68 + 0.1e1;
t73 = 0.1e1 / t76;
t67 = 0.1e1 / t71;
t63 = 0.1e1 / t82;
t62 = (0.1e1 + t90) * t88;
t59 = 0.1e1 / t61;
t58 = 0.1e1 / (t72 * t93 + 0.1e1);
t57 = (-t67 * t78 + t80 * t91) * t63 * t89;
t56 = (t76 * t59 - (-t79 * t92 + t65 * t75 + (-t65 * t87 + t92) * t62) * t75 * t60) * t81 * t58;
t1 = [t73 * t69 * t89, t62, 0, t62, t62, 0; (-t59 * t87 - (-t65 * t72 * t73 * t88 + (t69 - 0.1e1) * t75 * t64) * t75 * t93) * t58, t56, 0, t56, t56, 0; ((-t76 * t86 - t83) * t67 - (-t76 * t85 + t84) * t91) * t63, t57, 0, t57, t57, t82 * t63;];
Ja_rot  = t1;
