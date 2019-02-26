% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRP2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:25:21
% EndTime: 2019-02-26 21:25:21
% DurationCPUTime: 0.07s
% Computational Cost: add. (216->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
t51 = qJ(2) + pkin(9);
t49 = sin(t51);
t50 = cos(t51);
t53 = sin(qJ(1));
t61 = t53 * t50;
t44 = atan2(-t61, t49);
t42 = sin(t44);
t43 = cos(t44);
t35 = -t42 * t61 + t43 * t49;
t34 = 0.1e1 / t35 ^ 2;
t55 = cos(qJ(1));
t67 = t34 * t55 ^ 2;
t52 = sin(qJ(5));
t58 = t55 * t52;
t54 = cos(qJ(5));
t59 = t53 * t54;
t41 = t49 * t58 + t59;
t39 = 0.1e1 / t41 ^ 2;
t57 = t55 * t54;
t60 = t53 * t52;
t40 = -t49 * t57 + t60;
t66 = t39 * t40;
t65 = t42 * t49;
t48 = t50 ^ 2;
t64 = 0.1e1 / t49 ^ 2 * t48;
t63 = t50 * t55;
t45 = 0.1e1 / (t53 ^ 2 * t64 + 0.1e1);
t62 = t53 * t45;
t56 = t40 ^ 2 * t39 + 0.1e1;
t46 = 0.1e1 / t49;
t38 = 0.1e1 / t41;
t37 = 0.1e1 / t56;
t36 = (0.1e1 + t64) * t62;
t33 = 0.1e1 / t35;
t32 = 0.1e1 / (t48 * t67 + 0.1e1);
t1 = [-t46 * t45 * t63, t36, 0, 0, 0, 0; (-t33 * t61 - (t43 * t46 * t48 * t62 + (t45 - 0.1e1) * t50 * t42) * t50 * t67) * t32 (-t49 * t33 - (t53 * t65 + t43 * t50 + (-t43 * t61 - t65) * t36) * t50 * t34) * t55 * t32, 0, 0, 0, 0; ((t49 * t59 + t58) * t38 - (-t49 * t60 + t57) * t66) * t37 (-t38 * t54 - t52 * t66) * t37 * t63, 0, 0, t56 * t37, 0;];
Ja_rot  = t1;
