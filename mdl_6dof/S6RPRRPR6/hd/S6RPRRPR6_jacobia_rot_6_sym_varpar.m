% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:47
% EndTime: 2019-02-26 21:03:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (379->22), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->38)
t62 = pkin(10) + qJ(3);
t60 = cos(t62);
t59 = sin(t62);
t63 = sin(qJ(1));
t68 = t63 * t59;
t52 = atan2(-t68, -t60);
t50 = sin(t52);
t51 = cos(t52);
t44 = -t50 * t68 - t51 * t60;
t43 = 0.1e1 / t44 ^ 2;
t64 = cos(qJ(1));
t76 = t43 * t64 ^ 2;
t61 = qJ(4) + pkin(11) + qJ(6);
t55 = cos(t61);
t66 = t64 * t55;
t54 = sin(t61);
t70 = t63 * t54;
t49 = t60 * t66 + t70;
t47 = 0.1e1 / t49 ^ 2;
t67 = t64 * t54;
t69 = t63 * t55;
t48 = t60 * t67 - t69;
t75 = t47 * t48;
t74 = t50 * t60;
t56 = t59 ^ 2;
t73 = t56 / t60 ^ 2;
t72 = t59 * t64;
t53 = 0.1e1 / (t63 ^ 2 * t73 + 0.1e1);
t71 = t63 * t53;
t65 = t48 ^ 2 * t47 + 0.1e1;
t57 = 0.1e1 / t60;
t46 = 0.1e1 / t49;
t45 = (0.1e1 + t73) * t71;
t42 = 0.1e1 / t44;
t41 = 0.1e1 / t65;
t40 = 0.1e1 / (t56 * t76 + 0.1e1);
t39 = t65 * t41;
t1 = [t57 * t53 * t72, 0, t45, 0, 0, 0; (-t42 * t68 - (-t51 * t56 * t57 * t71 + (t53 - 0.1e1) * t59 * t50) * t59 * t76) * t40, 0 (t60 * t42 - (-t63 * t74 + t51 * t59 + (-t51 * t68 + t74) * t45) * t59 * t43) * t64 * t40, 0, 0, 0; ((-t60 * t70 - t66) * t46 - (-t60 * t69 + t67) * t75) * t41, 0 (-t46 * t54 + t55 * t75) * t41 * t72, t39, 0, t39;];
Ja_rot  = t1;
