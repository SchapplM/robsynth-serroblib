% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPPR4_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:00
% EndTime: 2019-02-26 20:00:00
% DurationCPUTime: 0.11s
% Computational Cost: add. (303->30), mult. (866->78), div. (60->9), fcn. (1237->13), ass. (0->48)
t61 = sin(pkin(10));
t64 = cos(pkin(10));
t69 = cos(qJ(2));
t65 = cos(pkin(6));
t67 = sin(qJ(2));
t72 = t65 * t67;
t54 = t61 * t69 + t64 * t72;
t66 = sin(qJ(3));
t62 = sin(pkin(6));
t68 = cos(qJ(3));
t74 = t62 * t68;
t46 = t54 * t66 + t64 * t74;
t75 = t62 * t66;
t57 = -t65 * t68 + t67 * t75;
t45 = atan2(-t46, t57);
t42 = sin(t45);
t43 = cos(t45);
t36 = -t42 * t46 + t43 * t57;
t35 = 0.1e1 / t36 ^ 2;
t56 = -t61 * t72 + t64 * t69;
t49 = t56 * t66 - t61 * t74;
t79 = t35 * t49;
t50 = t56 * t68 + t61 * t75;
t71 = t65 * t69;
t55 = t61 * t71 + t64 * t67;
t60 = sin(pkin(11));
t63 = cos(pkin(11));
t41 = t50 * t63 + t55 * t60;
t39 = 0.1e1 / t41 ^ 2;
t40 = t50 * t60 - t55 * t63;
t78 = t39 * t40;
t52 = 0.1e1 / t57 ^ 2;
t77 = t46 * t52;
t76 = t55 * t68;
t73 = t62 * t69;
t70 = -t42 * t57 - t43 * t46;
t58 = t65 * t66 + t67 * t74;
t53 = -t61 * t67 + t64 * t71;
t51 = 0.1e1 / t57;
t48 = t54 * t68 - t64 * t75;
t44 = 0.1e1 / (t46 ^ 2 * t52 + 0.1e1);
t38 = 0.1e1 / t41;
t37 = 0.1e1 / (t40 ^ 2 * t39 + 0.1e1);
t34 = 0.1e1 / t36;
t33 = 0.1e1 / (t49 ^ 2 * t35 + 0.1e1);
t32 = (-t51 * t53 + t73 * t77) * t66 * t44;
t31 = (-t48 * t51 + t58 * t77) * t44;
t1 = [0, t32, t31, 0, 0, 0; 0 (-t55 * t66 * t34 - ((-t42 * t53 + t43 * t73) * t66 + t70 * t32) * t79) * t33 (t50 * t34 - (t70 * t31 - t42 * t48 + t43 * t58) * t79) * t33, 0, 0, 0; 0 ((-t56 * t63 - t60 * t76) * t38 - (t56 * t60 - t63 * t76) * t78) * t37 (-t38 * t60 + t63 * t78) * t49 * t37, 0, 0, 0;];
Ja_rot  = t1;
