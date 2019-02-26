% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR2
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

function Ja_rot = S6RRPRRR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:54:40
% EndTime: 2019-02-26 21:54:40
% DurationCPUTime: 0.09s
% Computational Cost: add. (656->22), mult. (359->54), div. (84->9), fcn. (524->9), ass. (0->40)
t80 = qJ(2) + pkin(11) + qJ(4);
t79 = cos(t80);
t78 = sin(t80);
t84 = sin(qJ(1));
t91 = t84 * t78;
t73 = atan2(-t91, -t79);
t71 = sin(t73);
t72 = cos(t73);
t64 = -t71 * t91 - t72 * t79;
t63 = 0.1e1 / t64 ^ 2;
t85 = cos(qJ(1));
t97 = t63 * t85 ^ 2;
t83 = qJ(5) + qJ(6);
t82 = cos(t83);
t87 = t85 * t82;
t81 = sin(t83);
t90 = t84 * t81;
t70 = t79 * t87 + t90;
t68 = 0.1e1 / t70 ^ 2;
t88 = t85 * t81;
t89 = t84 * t82;
t69 = t79 * t88 - t89;
t96 = t68 * t69;
t95 = t71 * t79;
t75 = t78 ^ 2;
t94 = t75 / t79 ^ 2;
t93 = t78 * t85;
t74 = 0.1e1 / (t84 ^ 2 * t94 + 0.1e1);
t92 = t84 * t74;
t86 = t69 ^ 2 * t68 + 0.1e1;
t76 = 0.1e1 / t79;
t67 = 0.1e1 / t70;
t66 = 0.1e1 / t86;
t65 = (0.1e1 + t94) * t92;
t62 = 0.1e1 / t64;
t61 = 0.1e1 / (t75 * t97 + 0.1e1);
t60 = t86 * t66;
t59 = (-t67 * t81 + t82 * t96) * t66 * t93;
t58 = (t79 * t62 - (-t84 * t95 + t72 * t78 + (-t72 * t91 + t95) * t65) * t78 * t63) * t85 * t61;
t1 = [t76 * t74 * t93, t65, 0, t65, 0, 0; (-t62 * t91 - (-t72 * t75 * t76 * t92 + (t74 - 0.1e1) * t78 * t71) * t78 * t97) * t61, t58, 0, t58, 0, 0; ((-t79 * t90 - t87) * t67 - (-t79 * t89 + t88) * t96) * t66, t59, 0, t59, t60, t60;];
Ja_rot  = t1;
