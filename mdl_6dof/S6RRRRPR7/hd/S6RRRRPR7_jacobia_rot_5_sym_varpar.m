% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR7_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:06
% EndTime: 2019-02-26 22:34:06
% DurationCPUTime: 0.09s
% Computational Cost: add. (314->25), mult. (499->64), div. (77->11), fcn. (736->11), ass. (0->42)
t74 = sin(qJ(2));
t76 = cos(qJ(2));
t77 = cos(qJ(1));
t75 = sin(qJ(1));
t81 = cos(pkin(6));
t79 = t75 * t81;
t65 = -t74 * t79 + t77 * t76;
t69 = qJ(3) + qJ(4) + pkin(12);
t67 = sin(t69);
t68 = cos(t69);
t73 = sin(pkin(6));
t84 = t73 * t75;
t56 = t65 * t68 + t67 * t84;
t54 = 0.1e1 / t56 ^ 2;
t55 = t65 * t67 - t68 * t84;
t88 = t54 * t55;
t78 = t77 * t81;
t61 = t75 * t74 - t76 * t78;
t83 = t73 * t76;
t59 = atan2(-t61, -t83);
t58 = cos(t59);
t87 = t58 * t61;
t57 = sin(t59);
t51 = -t57 * t61 - t58 * t83;
t50 = 0.1e1 / t51 ^ 2;
t64 = t77 * t74 + t76 * t79;
t86 = t64 ^ 2 * t50;
t70 = 0.1e1 / t73;
t71 = 0.1e1 / t76;
t85 = t70 * t71;
t82 = t73 * t77;
t80 = t55 ^ 2 * t54 + 0.1e1;
t72 = 0.1e1 / t76 ^ 2;
t63 = t74 * t78 + t75 * t76;
t60 = 0.1e1 / (0.1e1 + t61 ^ 2 / t73 ^ 2 * t72);
t53 = 0.1e1 / t56;
t52 = 0.1e1 / t80;
t49 = 0.1e1 / t51;
t48 = 0.1e1 / (0.1e1 + t86);
t47 = (t61 * t72 * t74 + t63 * t71) * t70 * t60;
t46 = t80 * t52;
t1 = [t64 * t60 * t85, t47, 0, 0, 0, 0; (-t61 * t49 - (-t57 + (-t85 * t87 + t57) * t60) * t86) * t48 (t65 * t49 - (t58 * t73 * t74 - t57 * t63 + (t57 * t83 - t87) * t47) * t64 * t50) * t48, 0, 0, 0, 0; ((-t63 * t67 - t68 * t82) * t53 - (-t63 * t68 + t67 * t82) * t88) * t52 (-t67 * t53 + t68 * t88) * t64 * t52, t46, t46, 0, 0;];
Ja_rot  = t1;
