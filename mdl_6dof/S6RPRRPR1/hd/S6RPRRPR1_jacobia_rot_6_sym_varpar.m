% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:00:50
% EndTime: 2019-02-26 21:00:50
% DurationCPUTime: 0.09s
% Computational Cost: add. (713->22), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->39)
t77 = qJ(3) + qJ(4) + pkin(11);
t74 = cos(t77);
t73 = sin(t77);
t78 = qJ(1) + pkin(10);
t75 = sin(t78);
t86 = t75 * t73;
t68 = atan2(-t86, -t74);
t66 = sin(t68);
t67 = cos(t68);
t59 = -t66 * t86 - t67 * t74;
t58 = 0.1e1 / t59 ^ 2;
t76 = cos(t78);
t92 = t58 * t76 ^ 2;
t80 = cos(qJ(6));
t82 = t76 * t80;
t79 = sin(qJ(6));
t85 = t75 * t79;
t65 = t74 * t82 + t85;
t63 = 0.1e1 / t65 ^ 2;
t83 = t76 * t79;
t84 = t75 * t80;
t64 = t74 * t83 - t84;
t91 = t63 * t64;
t90 = t66 * t74;
t70 = t73 ^ 2;
t89 = t70 / t74 ^ 2;
t88 = t73 * t76;
t69 = 0.1e1 / (t75 ^ 2 * t89 + 0.1e1);
t87 = t75 * t69;
t81 = t64 ^ 2 * t63 + 0.1e1;
t71 = 0.1e1 / t74;
t62 = 0.1e1 / t65;
t61 = 0.1e1 / t81;
t60 = (0.1e1 + t89) * t87;
t57 = 0.1e1 / t59;
t56 = 0.1e1 / (t70 * t92 + 0.1e1);
t55 = (-t62 * t79 + t80 * t91) * t61 * t88;
t54 = (t74 * t57 - (-t75 * t90 + t67 * t73 + (-t67 * t86 + t90) * t60) * t73 * t58) * t76 * t56;
t1 = [t71 * t69 * t88, 0, t60, t60, 0, 0; (-t57 * t86 - (-t67 * t70 * t71 * t87 + (t69 - 0.1e1) * t73 * t66) * t73 * t92) * t56, 0, t54, t54, 0, 0; ((-t74 * t85 - t82) * t62 - (-t74 * t84 + t83) * t91) * t61, 0, t55, t55, 0, t81 * t61;];
Ja_rot  = t1;
