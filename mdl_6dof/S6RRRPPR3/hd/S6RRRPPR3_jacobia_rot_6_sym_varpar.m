% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:04:32
% EndTime: 2019-02-26 22:04:32
% DurationCPUTime: 0.11s
% Computational Cost: add. (324->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
t72 = qJ(2) + qJ(3);
t70 = sin(t72);
t71 = cos(t72);
t74 = sin(qJ(1));
t82 = t74 * t71;
t65 = atan2(-t82, t70);
t63 = sin(t65);
t64 = cos(t65);
t56 = -t63 * t82 + t64 * t70;
t55 = 0.1e1 / t56 ^ 2;
t76 = cos(qJ(1));
t88 = t55 * t76 ^ 2;
t75 = cos(qJ(6));
t78 = t76 * t75;
t73 = sin(qJ(6));
t81 = t74 * t73;
t62 = t70 * t78 - t81;
t60 = 0.1e1 / t62 ^ 2;
t79 = t76 * t73;
t80 = t74 * t75;
t61 = t70 * t79 + t80;
t87 = t60 * t61;
t86 = t63 * t70;
t69 = t71 ^ 2;
t85 = 0.1e1 / t70 ^ 2 * t69;
t84 = t71 * t76;
t66 = 0.1e1 / (t74 ^ 2 * t85 + 0.1e1);
t83 = t74 * t66;
t77 = t61 ^ 2 * t60 + 0.1e1;
t67 = 0.1e1 / t70;
t59 = 0.1e1 / t62;
t58 = 0.1e1 / t77;
t57 = (0.1e1 + t85) * t83;
t54 = 0.1e1 / t56;
t53 = 0.1e1 / (t69 * t88 + 0.1e1);
t52 = (t59 * t73 - t75 * t87) * t58 * t84;
t51 = (-t70 * t54 - (t74 * t86 + t64 * t71 + (-t64 * t82 - t86) * t57) * t71 * t55) * t76 * t53;
t1 = [-t67 * t66 * t84, t57, t57, 0, 0, 0; (-t54 * t82 - (t64 * t67 * t69 * t83 + (t66 - 0.1e1) * t71 * t63) * t71 * t88) * t53, t51, t51, 0, 0, 0; ((-t70 * t81 + t78) * t59 - (-t70 * t80 - t79) * t87) * t58, t52, t52, 0, 0, t77 * t58;];
Ja_rot  = t1;
