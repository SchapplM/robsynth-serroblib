% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:04:31
% EndTime: 2019-02-26 21:04:31
% DurationCPUTime: 0.11s
% Computational Cost: add. (519->20), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
t68 = qJ(3) + qJ(4) + pkin(10);
t66 = sin(t68);
t67 = cos(t68);
t72 = cos(qJ(1));
t76 = t72 * t67;
t61 = atan2(-t76, t66);
t55 = sin(t61);
t56 = cos(t61);
t52 = -t55 * t76 + t56 * t66;
t51 = 0.1e1 / t52 ^ 2;
t70 = sin(qJ(1));
t84 = t51 * t70 ^ 2;
t83 = t55 * t66;
t69 = sin(qJ(6));
t75 = t72 * t69;
t71 = cos(qJ(6));
t77 = t70 * t71;
t60 = t66 * t77 + t75;
t58 = 0.1e1 / t60 ^ 2;
t74 = t72 * t71;
t78 = t70 * t69;
t59 = t66 * t78 - t74;
t82 = t58 * t59;
t65 = t67 ^ 2;
t80 = 0.1e1 / t66 ^ 2 * t65;
t62 = 0.1e1 / (t72 ^ 2 * t80 + 0.1e1);
t81 = t62 * t72;
t79 = t67 * t70;
t73 = t59 ^ 2 * t58 + 0.1e1;
t63 = 0.1e1 / t66;
t57 = 0.1e1 / t60;
t54 = 0.1e1 / t73;
t53 = (0.1e1 + t80) * t81;
t50 = 0.1e1 / t52;
t49 = 0.1e1 / (t65 * t84 + 0.1e1);
t48 = (t57 * t69 - t71 * t82) * t54 * t79;
t47 = (t66 * t50 + (t72 * t83 + t56 * t67 + (-t56 * t76 - t83) * t53) * t67 * t51) * t70 * t49;
t1 = [t63 * t62 * t79, 0, t53, t53, 0, 0; (-t50 * t76 + (-t56 * t63 * t65 * t81 + (-t62 + 0.1e1) * t67 * t55) * t67 * t84) * t49, 0, t47, t47, 0, 0; ((t66 * t75 + t77) * t57 - (t66 * t74 - t78) * t82) * t54, 0, t48, t48, 0, t73 * t54;];
Ja_rot  = t1;
