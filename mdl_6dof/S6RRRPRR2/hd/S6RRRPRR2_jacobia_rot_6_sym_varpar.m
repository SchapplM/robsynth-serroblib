% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR2
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

function Ja_rot = S6RRRPRR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:16:33
% EndTime: 2019-02-26 22:16:33
% DurationCPUTime: 0.09s
% Computational Cost: add. (656->22), mult. (359->54), div. (84->9), fcn. (524->9), ass. (0->40)
t83 = qJ(2) + qJ(3) + pkin(11);
t82 = cos(t83);
t81 = sin(t83);
t87 = sin(qJ(1));
t94 = t87 * t81;
t76 = atan2(-t94, -t82);
t74 = sin(t76);
t75 = cos(t76);
t67 = -t74 * t94 - t75 * t82;
t66 = 0.1e1 / t67 ^ 2;
t88 = cos(qJ(1));
t100 = t66 * t88 ^ 2;
t86 = qJ(5) + qJ(6);
t85 = cos(t86);
t90 = t88 * t85;
t84 = sin(t86);
t93 = t87 * t84;
t73 = t82 * t90 + t93;
t71 = 0.1e1 / t73 ^ 2;
t91 = t88 * t84;
t92 = t87 * t85;
t72 = t82 * t91 - t92;
t99 = t71 * t72;
t98 = t74 * t82;
t78 = t81 ^ 2;
t97 = t78 / t82 ^ 2;
t96 = t81 * t88;
t77 = 0.1e1 / (t87 ^ 2 * t97 + 0.1e1);
t95 = t87 * t77;
t89 = t72 ^ 2 * t71 + 0.1e1;
t79 = 0.1e1 / t82;
t70 = 0.1e1 / t73;
t69 = 0.1e1 / t89;
t68 = (0.1e1 + t97) * t95;
t65 = 0.1e1 / t67;
t64 = 0.1e1 / (t78 * t100 + 0.1e1);
t63 = t89 * t69;
t62 = (-t70 * t84 + t85 * t99) * t69 * t96;
t61 = (t82 * t65 - (-t87 * t98 + t75 * t81 + (-t75 * t94 + t98) * t68) * t81 * t66) * t88 * t64;
t1 = [t79 * t77 * t96, t68, t68, 0, 0, 0; (-t65 * t94 - (-t75 * t78 * t79 * t95 + (t77 - 0.1e1) * t81 * t74) * t81 * t100) * t64, t61, t61, 0, 0, 0; ((-t82 * t93 - t90) * t70 - (-t82 * t92 + t91) * t99) * t69, t62, t62, 0, t63, t63;];
Ja_rot  = t1;
