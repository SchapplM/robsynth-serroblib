% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:48:18
% EndTime: 2019-02-26 22:48:18
% DurationCPUTime: 0.14s
% Computational Cost: add. (576->22), mult. (386->54), div. (89->9), fcn. (559->9), ass. (0->40)
t85 = qJ(2) + qJ(3);
t83 = cos(t85);
t82 = sin(t85);
t86 = sin(qJ(1));
t91 = t86 * t82;
t75 = atan2(-t91, -t83);
t73 = sin(t75);
t74 = cos(t75);
t67 = -t73 * t91 - t74 * t83;
t66 = 0.1e1 / t67 ^ 2;
t87 = cos(qJ(1));
t99 = t66 * t87 ^ 2;
t84 = qJ(4) + qJ(5) + qJ(6);
t78 = cos(t84);
t89 = t87 * t78;
t77 = sin(t84);
t93 = t86 * t77;
t72 = t83 * t89 + t93;
t70 = 0.1e1 / t72 ^ 2;
t90 = t87 * t77;
t92 = t86 * t78;
t71 = t83 * t90 - t92;
t98 = t70 * t71;
t97 = t73 * t83;
t79 = t82 ^ 2;
t96 = t79 / t83 ^ 2;
t95 = t82 * t87;
t76 = 0.1e1 / (t86 ^ 2 * t96 + 0.1e1);
t94 = t86 * t76;
t88 = t71 ^ 2 * t70 + 0.1e1;
t80 = 0.1e1 / t83;
t69 = 0.1e1 / t72;
t68 = (0.1e1 + t96) * t94;
t65 = 0.1e1 / t67;
t64 = 0.1e1 / t88;
t63 = 0.1e1 / (t79 * t99 + 0.1e1);
t62 = t88 * t64;
t61 = (-t69 * t77 + t78 * t98) * t64 * t95;
t60 = (t83 * t65 - (-t86 * t97 + t74 * t82 + (-t74 * t91 + t97) * t68) * t82 * t66) * t87 * t63;
t1 = [t80 * t76 * t95, t68, t68, 0, 0, 0; (-t65 * t91 - (-t74 * t79 * t80 * t94 + (t76 - 0.1e1) * t82 * t73) * t82 * t99) * t63, t60, t60, 0, 0, 0; ((-t83 * t93 - t89) * t69 - (-t83 * t92 + t90) * t98) * t64, t61, t61, t62, t62, t62;];
Ja_rot  = t1;
