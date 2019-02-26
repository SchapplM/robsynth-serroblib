% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:16:03
% EndTime: 2019-02-26 22:16:03
% DurationCPUTime: 0.10s
% Computational Cost: add. (1004->21), mult. (440->54), div. (106->9), fcn. (646->9), ass. (0->38)
t78 = qJ(2) + qJ(3) + pkin(11) + qJ(5);
t77 = cos(t78);
t76 = sin(t78);
t80 = sin(qJ(1));
t88 = t80 * t76;
t67 = atan2(-t88, -t77);
t65 = sin(t67);
t66 = cos(t67);
t62 = -t65 * t88 - t66 * t77;
t61 = 0.1e1 / t62 ^ 2;
t82 = cos(qJ(1));
t94 = t61 * t82 ^ 2;
t93 = t65 * t77;
t81 = cos(qJ(6));
t84 = t82 * t81;
t79 = sin(qJ(6));
t87 = t80 * t79;
t72 = t77 * t84 + t87;
t69 = 0.1e1 / t72 ^ 2;
t85 = t82 * t79;
t86 = t80 * t81;
t71 = t77 * t85 - t86;
t92 = t69 * t71;
t73 = t76 ^ 2;
t91 = t73 / t77 ^ 2;
t90 = t76 * t82;
t70 = 0.1e1 / (t80 ^ 2 * t91 + 0.1e1);
t89 = t80 * t70;
t83 = t71 ^ 2 * t69 + 0.1e1;
t74 = 0.1e1 / t77;
t68 = 0.1e1 / t72;
t64 = 0.1e1 / t83;
t63 = (0.1e1 + t91) * t89;
t60 = 0.1e1 / t62;
t59 = 0.1e1 / (t73 * t94 + 0.1e1);
t58 = (-t68 * t79 + t81 * t92) * t64 * t90;
t57 = (t77 * t60 - (-t80 * t93 + t66 * t76 + (-t66 * t88 + t93) * t63) * t76 * t61) * t82 * t59;
t1 = [t74 * t70 * t90, t63, t63, 0, t63, 0; (-t60 * t88 - (-t66 * t73 * t74 * t89 + (t70 - 0.1e1) * t76 * t65) * t76 * t94) * t59, t57, t57, 0, t57, 0; ((-t77 * t87 - t84) * t68 - (-t77 * t86 + t85) * t92) * t64, t58, t58, 0, t58, t83 * t64;];
Ja_rot  = t1;
