% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR3
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
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:10
% EndTime: 2019-02-26 21:55:10
% DurationCPUTime: 0.08s
% Computational Cost: add. (424->22), mult. (278->54), div. (62->9), fcn. (402->9), ass. (0->38)
t71 = qJ(2) + pkin(11);
t69 = cos(t71);
t68 = sin(t71);
t72 = sin(qJ(1));
t77 = t72 * t68;
t61 = atan2(-t77, -t69);
t59 = sin(t61);
t60 = cos(t61);
t53 = -t59 * t77 - t60 * t69;
t51 = 0.1e1 / t53 ^ 2;
t73 = cos(qJ(1));
t85 = t51 * t73 ^ 2;
t70 = qJ(4) + qJ(5) + qJ(6);
t64 = cos(t70);
t75 = t73 * t64;
t63 = sin(t70);
t79 = t72 * t63;
t58 = t69 * t75 + t79;
t56 = 0.1e1 / t58 ^ 2;
t76 = t73 * t63;
t78 = t72 * t64;
t57 = t69 * t76 - t78;
t84 = t56 * t57;
t83 = t59 * t69;
t65 = t68 ^ 2;
t82 = t65 / t69 ^ 2;
t81 = t68 * t73;
t62 = 0.1e1 / (t72 ^ 2 * t82 + 0.1e1);
t80 = t72 * t62;
t74 = t57 ^ 2 * t56 + 0.1e1;
t66 = 0.1e1 / t69;
t55 = 0.1e1 / t58;
t54 = (0.1e1 + t82) * t80;
t52 = 0.1e1 / t74;
t50 = 0.1e1 / t53;
t49 = 0.1e1 / (t65 * t85 + 0.1e1);
t48 = t74 * t52;
t1 = [t66 * t62 * t81, t54, 0, 0, 0, 0; (-t50 * t77 - (-t60 * t65 * t66 * t80 + (t62 - 0.1e1) * t68 * t59) * t68 * t85) * t49 (t69 * t50 - (-t72 * t83 + t60 * t68 + (-t60 * t77 + t83) * t54) * t68 * t51) * t73 * t49, 0, 0, 0, 0; ((-t69 * t79 - t75) * t55 - (-t69 * t78 + t76) * t84) * t52 (-t55 * t63 + t64 * t84) * t52 * t81, 0, t48, t48, t48;];
Ja_rot  = t1;
