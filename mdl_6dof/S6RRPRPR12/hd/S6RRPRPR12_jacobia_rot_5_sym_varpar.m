% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR12_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:12
% EndTime: 2019-02-26 21:44:12
% DurationCPUTime: 0.09s
% Computational Cost: add. (197->25), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->41)
t62 = sin(qJ(2));
t64 = cos(qJ(2));
t65 = cos(qJ(1));
t63 = sin(qJ(1));
t69 = cos(pkin(6));
t67 = t63 * t69;
t51 = t65 * t62 + t64 * t67;
t58 = qJ(4) + pkin(11);
t55 = sin(t58);
t56 = cos(t58);
t61 = sin(pkin(6));
t71 = t61 * t63;
t43 = t51 * t55 + t56 * t71;
t41 = 0.1e1 / t43 ^ 2;
t42 = -t51 * t56 + t55 * t71;
t76 = t41 * t42;
t66 = t65 * t69;
t49 = t62 * t66 + t63 * t64;
t72 = t61 * t62;
t47 = atan2(-t49, t72);
t45 = cos(t47);
t75 = t45 * t49;
t44 = sin(t47);
t38 = -t44 * t49 + t45 * t72;
t37 = 0.1e1 / t38 ^ 2;
t52 = -t62 * t67 + t65 * t64;
t74 = t52 ^ 2 * t37;
t57 = 0.1e1 / t61;
t59 = 0.1e1 / t62;
t73 = t57 * t59;
t70 = t61 * t65;
t68 = t42 ^ 2 * t41 + 0.1e1;
t60 = 0.1e1 / t62 ^ 2;
t48 = t63 * t62 - t64 * t66;
t46 = 0.1e1 / (0.1e1 + t49 ^ 2 / t61 ^ 2 * t60);
t40 = 0.1e1 / t43;
t39 = 0.1e1 / t68;
t36 = 0.1e1 / t38;
t35 = 0.1e1 / (0.1e1 + t74);
t34 = (t49 * t60 * t64 + t48 * t59) * t57 * t46;
t1 = [-t52 * t46 * t73, t34, 0, 0, 0, 0; (-t49 * t36 - (-t44 + (t73 * t75 + t44) * t46) * t74) * t35 (-t51 * t36 - (t45 * t61 * t64 + t44 * t48 + (-t44 * t72 - t75) * t34) * t52 * t37) * t35, 0, 0, 0, 0; ((t48 * t56 + t55 * t70) * t40 - (-t48 * t55 + t56 * t70) * t76) * t39 (-t56 * t40 - t55 * t76) * t52 * t39, 0, t68 * t39, 0, 0;];
Ja_rot  = t1;
