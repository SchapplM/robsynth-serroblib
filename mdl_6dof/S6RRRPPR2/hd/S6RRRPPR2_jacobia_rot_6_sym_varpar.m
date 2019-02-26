% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:56
% EndTime: 2019-02-26 22:03:56
% DurationCPUTime: 0.12s
% Computational Cost: add. (520->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
t76 = qJ(2) + qJ(3) + pkin(10);
t74 = sin(t76);
t75 = cos(t76);
t78 = sin(qJ(1));
t86 = t78 * t75;
t69 = atan2(-t86, t74);
t63 = sin(t69);
t64 = cos(t69);
t60 = -t63 * t86 + t64 * t74;
t59 = 0.1e1 / t60 ^ 2;
t80 = cos(qJ(1));
t92 = t59 * t80 ^ 2;
t91 = t63 * t74;
t77 = sin(qJ(6));
t83 = t80 * t77;
t79 = cos(qJ(6));
t84 = t78 * t79;
t68 = t74 * t83 + t84;
t66 = 0.1e1 / t68 ^ 2;
t82 = t80 * t79;
t85 = t78 * t77;
t67 = -t74 * t82 + t85;
t90 = t66 * t67;
t73 = t75 ^ 2;
t89 = 0.1e1 / t74 ^ 2 * t73;
t88 = t75 * t80;
t70 = 0.1e1 / (t78 ^ 2 * t89 + 0.1e1);
t87 = t78 * t70;
t81 = t67 ^ 2 * t66 + 0.1e1;
t71 = 0.1e1 / t74;
t65 = 0.1e1 / t68;
t62 = 0.1e1 / t81;
t61 = (0.1e1 + t89) * t87;
t58 = 0.1e1 / t60;
t57 = 0.1e1 / (t73 * t92 + 0.1e1);
t56 = (-t65 * t79 - t77 * t90) * t62 * t88;
t55 = (-t74 * t58 - (t78 * t91 + t64 * t75 + (-t64 * t86 - t91) * t61) * t75 * t59) * t80 * t57;
t1 = [-t71 * t70 * t88, t61, t61, 0, 0, 0; (-t58 * t86 - (t64 * t71 * t73 * t87 + (t70 - 0.1e1) * t75 * t63) * t75 * t92) * t57, t55, t55, 0, 0, 0; ((t74 * t84 + t83) * t65 - (-t74 * t85 + t82) * t90) * t62, t56, t56, 0, 0, t81 * t62;];
Ja_rot  = t1;
