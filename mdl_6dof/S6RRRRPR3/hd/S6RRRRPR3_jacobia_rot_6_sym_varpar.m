% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:42
% EndTime: 2019-02-26 22:31:42
% DurationCPUTime: 0.09s
% Computational Cost: add. (695->21), mult. (440->54), div. (106->9), fcn. (646->9), ass. (0->38)
t77 = qJ(2) + qJ(3) + qJ(4);
t75 = sin(t77);
t76 = cos(t77);
t79 = sin(qJ(1));
t87 = t79 * t76;
t70 = atan2(-t87, t75);
t66 = sin(t70);
t67 = cos(t70);
t61 = -t66 * t87 + t67 * t75;
t60 = 0.1e1 / t61 ^ 2;
t81 = cos(qJ(1));
t93 = t60 * t81 ^ 2;
t78 = sin(qJ(6));
t84 = t81 * t78;
t80 = cos(qJ(6));
t85 = t79 * t80;
t69 = t75 * t84 + t85;
t65 = 0.1e1 / t69 ^ 2;
t83 = t81 * t80;
t86 = t79 * t78;
t68 = -t75 * t83 + t86;
t92 = t65 * t68;
t91 = t66 * t75;
t74 = t76 ^ 2;
t90 = 0.1e1 / t75 ^ 2 * t74;
t89 = t76 * t81;
t71 = 0.1e1 / (t79 ^ 2 * t90 + 0.1e1);
t88 = t79 * t71;
t82 = t68 ^ 2 * t65 + 0.1e1;
t72 = 0.1e1 / t75;
t64 = 0.1e1 / t69;
t63 = 0.1e1 / t82;
t62 = (0.1e1 + t90) * t88;
t59 = 0.1e1 / t61;
t58 = 0.1e1 / (t74 * t93 + 0.1e1);
t57 = (-t64 * t80 - t78 * t92) * t63 * t89;
t56 = (-t75 * t59 - (t79 * t91 + t67 * t76 + (-t67 * t87 - t91) * t62) * t76 * t60) * t81 * t58;
t1 = [-t72 * t71 * t89, t62, t62, t62, 0, 0; (-t59 * t87 - (t67 * t72 * t74 * t88 + (t71 - 0.1e1) * t76 * t66) * t76 * t93) * t58, t56, t56, t56, 0, 0; ((t75 * t85 + t84) * t64 - (-t75 * t86 + t83) * t92) * t63, t57, t57, t57, 0, t82 * t63;];
Ja_rot  = t1;
