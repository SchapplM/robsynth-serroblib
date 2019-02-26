% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRP1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:50
% EndTime: 2019-02-26 21:24:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (640->28), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->42)
t65 = qJ(2) + pkin(9);
t61 = sin(t65);
t81 = t61 ^ 2;
t63 = cos(t65);
t64 = pkin(10) + qJ(5);
t62 = cos(t64);
t68 = cos(qJ(1));
t70 = t68 * t62;
t60 = sin(t64);
t67 = sin(qJ(1));
t73 = t67 * t60;
t48 = t63 * t73 + t70;
t75 = t61 * t60;
t44 = atan2(-t48, t75);
t41 = sin(t44);
t42 = cos(t44);
t40 = -t41 * t48 + t42 * t75;
t39 = 0.1e1 / t40 ^ 2;
t71 = t68 * t60;
t72 = t67 * t62;
t51 = t63 * t71 - t72;
t80 = t39 * t51;
t78 = t42 * t48;
t77 = t51 ^ 2 * t39;
t55 = 0.1e1 / t60;
t58 = 0.1e1 / t61;
t76 = t55 * t58;
t74 = t61 * t68;
t52 = t63 * t70 + t73;
t47 = 0.1e1 / t52 ^ 2;
t69 = t68 ^ 2 * t81 * t47;
t59 = 0.1e1 / t81;
t56 = 0.1e1 / t60 ^ 2;
t50 = t63 * t72 - t71;
t46 = 0.1e1 / t52;
t45 = 0.1e1 / (0.1e1 + t69);
t43 = 0.1e1 / (t48 ^ 2 * t59 * t56 + 0.1e1);
t38 = 0.1e1 / t40;
t37 = (t48 * t55 * t59 * t63 + t67) * t43;
t36 = 0.1e1 / (0.1e1 + t77);
t35 = (t48 * t56 * t62 - t50 * t55) * t58 * t43;
t1 = [-t51 * t43 * t76, t37, 0, 0, t35, 0; (-t48 * t38 - (-t41 + (t76 * t78 + t41) * t43) * t77) * t36 (t37 * t78 * t80 + (-t38 * t74 - (t42 * t63 + (-t37 + t67) * t61 * t41) * t80) * t60) * t36, 0, 0 (t52 * t38 - (t42 * t61 * t62 - t41 * t50 + (-t41 * t75 - t78) * t35) * t80) * t36, 0; (-t47 * t50 * t68 + t46 * t67) * t61 * t45 (-t46 * t63 * t68 - t62 * t69) * t45, 0, 0, -t51 * t47 * t45 * t74, 0;];
Ja_rot  = t1;
