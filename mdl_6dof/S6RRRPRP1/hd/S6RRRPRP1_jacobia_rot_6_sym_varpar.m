% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:09:30
% EndTime: 2019-02-26 22:09:30
% DurationCPUTime: 0.08s
% Computational Cost: add. (554->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
t72 = qJ(2) + qJ(3) + pkin(10);
t71 = cos(t72);
t70 = sin(t72);
t74 = sin(qJ(1));
t82 = t74 * t70;
t63 = atan2(-t82, -t71);
t59 = sin(t63);
t60 = cos(t63);
t56 = -t59 * t82 - t60 * t71;
t55 = 0.1e1 / t56 ^ 2;
t76 = cos(qJ(1));
t88 = t55 * t76 ^ 2;
t87 = t59 * t71;
t75 = cos(qJ(5));
t78 = t76 * t75;
t73 = sin(qJ(5));
t81 = t74 * t73;
t65 = t71 * t78 + t81;
t62 = 0.1e1 / t65 ^ 2;
t79 = t76 * t73;
t80 = t74 * t75;
t64 = t71 * t79 - t80;
t86 = t62 * t64;
t67 = t70 ^ 2;
t85 = t67 / t71 ^ 2;
t84 = t70 * t76;
t66 = 0.1e1 / (t74 ^ 2 * t85 + 0.1e1);
t83 = t74 * t66;
t77 = t64 ^ 2 * t62 + 0.1e1;
t68 = 0.1e1 / t71;
t61 = 0.1e1 / t65;
t58 = 0.1e1 / t77;
t57 = (0.1e1 + t85) * t83;
t54 = 0.1e1 / t56;
t53 = 0.1e1 / (t67 * t88 + 0.1e1);
t52 = (-t61 * t73 + t75 * t86) * t58 * t84;
t51 = (t71 * t54 - (-t74 * t87 + t60 * t70 + (-t60 * t82 + t87) * t57) * t70 * t55) * t76 * t53;
t1 = [t68 * t66 * t84, t57, t57, 0, 0, 0; (-t54 * t82 - (-t60 * t67 * t68 * t83 + (t66 - 0.1e1) * t70 * t59) * t70 * t88) * t53, t51, t51, 0, 0, 0; ((-t71 * t81 - t78) * t61 - (-t71 * t80 + t79) * t86) * t58, t52, t52, 0, t77 * t58, 0;];
Ja_rot  = t1;
