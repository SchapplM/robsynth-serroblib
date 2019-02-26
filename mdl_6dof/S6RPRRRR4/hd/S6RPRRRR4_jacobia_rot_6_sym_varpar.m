% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:16:26
% EndTime: 2019-02-26 21:16:26
% DurationCPUTime: 0.09s
% Computational Cost: add. (1004->21), mult. (440->54), div. (106->9), fcn. (646->9), ass. (0->38)
t75 = pkin(11) + qJ(3) + qJ(4) + qJ(5);
t74 = cos(t75);
t73 = sin(t75);
t77 = sin(qJ(1));
t85 = t77 * t73;
t64 = atan2(-t85, -t74);
t62 = sin(t64);
t63 = cos(t64);
t59 = -t62 * t85 - t63 * t74;
t58 = 0.1e1 / t59 ^ 2;
t79 = cos(qJ(1));
t91 = t58 * t79 ^ 2;
t90 = t62 * t74;
t78 = cos(qJ(6));
t81 = t79 * t78;
t76 = sin(qJ(6));
t84 = t77 * t76;
t69 = t74 * t81 + t84;
t66 = 0.1e1 / t69 ^ 2;
t82 = t79 * t76;
t83 = t77 * t78;
t68 = t74 * t82 - t83;
t89 = t66 * t68;
t70 = t73 ^ 2;
t88 = t70 / t74 ^ 2;
t87 = t73 * t79;
t67 = 0.1e1 / (t77 ^ 2 * t88 + 0.1e1);
t86 = t77 * t67;
t80 = t68 ^ 2 * t66 + 0.1e1;
t71 = 0.1e1 / t74;
t65 = 0.1e1 / t69;
t61 = 0.1e1 / t80;
t60 = (0.1e1 + t88) * t86;
t57 = 0.1e1 / t59;
t56 = 0.1e1 / (t70 * t91 + 0.1e1);
t55 = (-t65 * t76 + t78 * t89) * t61 * t87;
t54 = (t74 * t57 - (-t77 * t90 + t63 * t73 + (-t63 * t85 + t90) * t60) * t73 * t58) * t79 * t56;
t1 = [t71 * t67 * t87, 0, t60, t60, t60, 0; (-t57 * t85 - (-t63 * t70 * t71 * t86 + (t67 - 0.1e1) * t73 * t62) * t73 * t91) * t56, 0, t54, t54, t54, 0; ((-t74 * t84 - t81) * t65 - (-t74 * t83 + t82) * t89) * t61, 0, t55, t55, t55, t80 * t61;];
Ja_rot  = t1;
