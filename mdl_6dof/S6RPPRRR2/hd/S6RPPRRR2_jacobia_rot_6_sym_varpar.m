% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:35:24
% EndTime: 2019-02-26 20:35:25
% DurationCPUTime: 0.10s
% Computational Cost: add. (440->23), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->39)
t64 = pkin(11) + qJ(4);
t60 = cos(t64);
t58 = sin(t64);
t65 = qJ(1) + pkin(10);
t59 = sin(t65);
t72 = t59 * t58;
t53 = atan2(-t72, -t60);
t51 = sin(t53);
t52 = cos(t53);
t44 = -t51 * t72 - t52 * t60;
t43 = 0.1e1 / t44 ^ 2;
t61 = cos(t65);
t78 = t43 * t61 ^ 2;
t66 = qJ(5) + qJ(6);
t63 = cos(t66);
t68 = t61 * t63;
t62 = sin(t66);
t71 = t59 * t62;
t50 = t60 * t68 + t71;
t48 = 0.1e1 / t50 ^ 2;
t69 = t61 * t62;
t70 = t59 * t63;
t49 = t60 * t69 - t70;
t77 = t48 * t49;
t76 = t51 * t60;
t55 = t58 ^ 2;
t75 = t55 / t60 ^ 2;
t74 = t58 * t61;
t54 = 0.1e1 / (t59 ^ 2 * t75 + 0.1e1);
t73 = t59 * t54;
t67 = t49 ^ 2 * t48 + 0.1e1;
t56 = 0.1e1 / t60;
t47 = 0.1e1 / t50;
t46 = 0.1e1 / t67;
t45 = (0.1e1 + t75) * t73;
t42 = 0.1e1 / t44;
t41 = 0.1e1 / (t55 * t78 + 0.1e1);
t40 = t67 * t46;
t1 = [t56 * t54 * t74, 0, 0, t45, 0, 0; (-t42 * t72 - (-t52 * t55 * t56 * t73 + (t54 - 0.1e1) * t58 * t51) * t58 * t78) * t41, 0, 0 (t60 * t42 - (-t59 * t76 + t52 * t58 + (-t52 * t72 + t76) * t45) * t58 * t43) * t61 * t41, 0, 0; ((-t60 * t71 - t68) * t47 - (-t60 * t70 + t69) * t77) * t46, 0, 0 (-t47 * t62 + t63 * t77) * t46 * t74, t40, t40;];
Ja_rot  = t1;
