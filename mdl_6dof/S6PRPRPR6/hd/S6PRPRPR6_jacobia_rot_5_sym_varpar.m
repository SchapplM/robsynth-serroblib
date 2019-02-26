% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR6_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:12
% EndTime: 2019-02-26 19:49:12
% DurationCPUTime: 0.11s
% Computational Cost: add. (303->29), mult. (866->79), div. (60->9), fcn. (1237->13), ass. (0->49)
t61 = sin(pkin(10));
t64 = cos(pkin(10));
t67 = sin(qJ(2));
t65 = cos(pkin(6));
t69 = cos(qJ(2));
t71 = t65 * t69;
t54 = t61 * t67 - t64 * t71;
t68 = cos(qJ(4));
t62 = sin(pkin(6));
t66 = sin(qJ(4));
t76 = t62 * t66;
t49 = t54 * t68 + t64 * t76;
t73 = t62 * t69;
t58 = t65 * t66 + t68 * t73;
t46 = atan2(t49, t58);
t43 = sin(t46);
t44 = cos(t46);
t37 = t43 * t49 + t44 * t58;
t36 = 0.1e1 / t37 ^ 2;
t56 = t61 * t71 + t64 * t67;
t47 = -t56 * t68 + t61 * t76;
t80 = t36 * t47;
t74 = t62 * t68;
t48 = t56 * t66 + t61 * t74;
t72 = t65 * t67;
t57 = -t61 * t72 + t64 * t69;
t60 = sin(pkin(11));
t63 = cos(pkin(11));
t42 = t48 * t63 + t57 * t60;
t40 = 0.1e1 / t42 ^ 2;
t41 = t48 * t60 - t57 * t63;
t79 = t40 * t41;
t53 = 0.1e1 / t58 ^ 2;
t78 = t49 * t53;
t77 = t57 * t66;
t75 = t62 * t67;
t70 = -t43 * t58 + t44 * t49;
t59 = t65 * t68 - t66 * t73;
t55 = t61 * t69 + t64 * t72;
t52 = 0.1e1 / t58;
t50 = -t54 * t66 + t64 * t74;
t45 = 0.1e1 / (t49 ^ 2 * t53 + 0.1e1);
t39 = 0.1e1 / t42;
t38 = 0.1e1 / (t41 ^ 2 * t40 + 0.1e1);
t35 = 0.1e1 / t37;
t34 = 0.1e1 / (t47 ^ 2 * t36 + 0.1e1);
t33 = (t52 * t55 + t75 * t78) * t68 * t45;
t32 = (t50 * t52 - t59 * t78) * t45;
t1 = [0, t33, 0, t32, 0, 0; 0 (-t57 * t68 * t35 - ((t43 * t55 - t44 * t75) * t68 + t70 * t33) * t80) * t34, 0 (t48 * t35 - (t70 * t32 + t43 * t50 + t44 * t59) * t80) * t34, 0, 0; 0 ((t56 * t63 + t60 * t77) * t39 - (-t56 * t60 + t63 * t77) * t79) * t38, 0 (-t39 * t60 + t63 * t79) * t47 * t38, 0, 0;];
Ja_rot  = t1;
