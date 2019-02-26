% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:40
% EndTime: 2019-02-26 21:37:40
% DurationCPUTime: 0.09s
% Computational Cost: add. (530->21), mult. (305->53), div. (74->9), fcn. (454->9), ass. (0->37)
t65 = qJ(2) + pkin(10) + qJ(4);
t64 = cos(t65);
t63 = sin(t65);
t68 = sin(qJ(1));
t74 = t68 * t63;
t58 = atan2(-t74, -t64);
t52 = sin(t58);
t53 = cos(t58);
t49 = -t52 * t74 - t53 * t64;
t48 = 0.1e1 / t49 ^ 2;
t69 = cos(qJ(1));
t80 = t48 * t69 ^ 2;
t79 = t52 * t64;
t67 = cos(pkin(11));
t70 = t69 * t67;
t66 = sin(pkin(11));
t73 = t68 * t66;
t57 = t64 * t70 + t73;
t55 = 0.1e1 / t57 ^ 2;
t71 = t69 * t66;
t72 = t68 * t67;
t56 = t64 * t71 - t72;
t78 = t55 * t56;
t60 = t63 ^ 2;
t77 = t60 / t64 ^ 2;
t76 = t63 * t69;
t59 = 0.1e1 / (t68 ^ 2 * t77 + 0.1e1);
t75 = t68 * t59;
t61 = 0.1e1 / t64;
t54 = 0.1e1 / t57;
t51 = 0.1e1 / (t56 ^ 2 * t55 + 0.1e1);
t50 = (0.1e1 + t77) * t75;
t47 = 0.1e1 / t49;
t46 = 0.1e1 / (t60 * t80 + 0.1e1);
t45 = (-t54 * t66 + t67 * t78) * t51 * t76;
t44 = (t64 * t47 - (-t68 * t79 + t53 * t63 + (-t53 * t74 + t79) * t50) * t63 * t48) * t69 * t46;
t1 = [t61 * t59 * t76, t50, 0, t50, 0, 0; (-t47 * t74 - (-t53 * t60 * t61 * t75 + (t59 - 0.1e1) * t63 * t52) * t63 * t80) * t46, t44, 0, t44, 0, 0; ((-t64 * t73 - t70) * t54 - (-t64 * t72 + t71) * t78) * t51, t45, 0, t45, 0, 0;];
Ja_rot  = t1;
