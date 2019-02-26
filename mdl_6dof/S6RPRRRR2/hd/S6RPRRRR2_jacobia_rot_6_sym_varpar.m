% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR2
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
% Datum: 2019-02-26 21:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:15:25
% EndTime: 2019-02-26 21:15:25
% DurationCPUTime: 0.09s
% Computational Cost: add. (626->23), mult. (359->57), div. (84->9), fcn. (524->9), ass. (0->40)
t79 = qJ(1) + pkin(11);
t74 = cos(t79);
t92 = t74 ^ 2;
t81 = qJ(3) + qJ(4);
t78 = cos(t81);
t73 = sin(t79);
t76 = sin(t81);
t86 = t73 * t76;
t68 = atan2(-t86, -t78);
t66 = sin(t68);
t67 = cos(t68);
t59 = -t66 * t86 - t67 * t78;
t58 = 0.1e1 / t59 ^ 2;
t91 = t58 * t76;
t80 = qJ(5) + qJ(6);
t75 = sin(t80);
t77 = cos(t80);
t83 = t77 * t78;
t65 = t73 * t75 + t74 * t83;
t63 = 0.1e1 / t65 ^ 2;
t84 = t75 * t78;
t64 = -t73 * t77 + t74 * t84;
t90 = t63 * t64;
t89 = t66 * t78;
t70 = t76 ^ 2;
t88 = t70 / t78 ^ 2;
t69 = 0.1e1 / (t73 ^ 2 * t88 + 0.1e1);
t87 = t73 * t69;
t85 = t74 * t76;
t82 = t64 ^ 2 * t63 + 0.1e1;
t71 = 0.1e1 / t78;
t62 = 0.1e1 / t65;
t61 = (0.1e1 + t88) * t87;
t60 = 0.1e1 / t82;
t57 = 0.1e1 / t59;
t56 = 0.1e1 / (t92 * t70 * t58 + 0.1e1);
t55 = t82 * t60;
t54 = (-t62 * t75 + t77 * t90) * t60 * t85;
t53 = (t78 * t57 - (-t73 * t89 + t67 * t76 + (-t67 * t86 + t89) * t61) * t91) * t74 * t56;
t1 = [t71 * t69 * t85, 0, t61, t61, 0, 0; (-t57 * t86 - (-t67 * t70 * t71 * t87 + (t69 - 0.1e1) * t76 * t66) * t92 * t91) * t56, 0, t53, t53, 0, 0; ((-t73 * t84 - t74 * t77) * t62 - (-t73 * t83 + t74 * t75) * t90) * t60, 0, t54, t54, t55, t55;];
Ja_rot  = t1;
