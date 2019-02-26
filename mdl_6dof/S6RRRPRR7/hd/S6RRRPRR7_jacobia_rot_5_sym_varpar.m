% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR7
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR7_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:18
% EndTime: 2019-02-26 22:19:18
% DurationCPUTime: 0.12s
% Computational Cost: add. (314->25), mult. (499->64), div. (77->11), fcn. (736->11), ass. (0->42)
t72 = sin(qJ(2));
t74 = cos(qJ(2));
t75 = cos(qJ(1));
t73 = sin(qJ(1));
t79 = cos(pkin(6));
t77 = t73 * t79;
t63 = -t72 * t77 + t75 * t74;
t67 = qJ(3) + pkin(12) + qJ(5);
t65 = sin(t67);
t66 = cos(t67);
t71 = sin(pkin(6));
t82 = t71 * t73;
t54 = t63 * t66 + t65 * t82;
t52 = 0.1e1 / t54 ^ 2;
t53 = t63 * t65 - t66 * t82;
t86 = t52 * t53;
t76 = t75 * t79;
t59 = t73 * t72 - t74 * t76;
t81 = t71 * t74;
t57 = atan2(-t59, -t81);
t56 = cos(t57);
t85 = t56 * t59;
t55 = sin(t57);
t49 = -t55 * t59 - t56 * t81;
t48 = 0.1e1 / t49 ^ 2;
t62 = t75 * t72 + t74 * t77;
t84 = t62 ^ 2 * t48;
t68 = 0.1e1 / t71;
t69 = 0.1e1 / t74;
t83 = t68 * t69;
t80 = t71 * t75;
t78 = t53 ^ 2 * t52 + 0.1e1;
t70 = 0.1e1 / t74 ^ 2;
t61 = t72 * t76 + t73 * t74;
t58 = 0.1e1 / (0.1e1 + t59 ^ 2 / t71 ^ 2 * t70);
t51 = 0.1e1 / t54;
t50 = 0.1e1 / t78;
t47 = 0.1e1 / t49;
t46 = 0.1e1 / (0.1e1 + t84);
t45 = (t59 * t70 * t72 + t61 * t69) * t68 * t58;
t44 = t78 * t50;
t1 = [t62 * t58 * t83, t45, 0, 0, 0, 0; (-t59 * t47 - (-t55 + (-t83 * t85 + t55) * t58) * t84) * t46 (t63 * t47 - (t56 * t71 * t72 - t55 * t61 + (t55 * t81 - t85) * t45) * t62 * t48) * t46, 0, 0, 0, 0; ((-t61 * t65 - t66 * t80) * t51 - (-t61 * t66 + t65 * t80) * t86) * t50 (-t65 * t51 + t66 * t86) * t62 * t50, t44, 0, t44, 0;];
Ja_rot  = t1;
