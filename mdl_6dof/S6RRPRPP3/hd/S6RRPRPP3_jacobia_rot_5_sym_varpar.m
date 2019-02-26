% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPP3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:55
% EndTime: 2019-02-26 21:35:56
% DurationCPUTime: 0.10s
% Computational Cost: add. (437->27), mult. (486->70), div. (108->11), fcn. (755->9), ass. (0->39)
t50 = pkin(9) + qJ(4);
t48 = sin(t50);
t49 = cos(t50);
t58 = cos(qJ(1));
t59 = t58 * t49;
t56 = sin(qJ(1));
t57 = cos(qJ(2));
t61 = t56 * t57;
t39 = t48 * t61 + t59;
t55 = sin(qJ(2));
t62 = t55 * t48;
t37 = atan2(-t39, t62);
t34 = sin(t37);
t35 = cos(t37);
t33 = -t34 * t39 + t35 * t62;
t32 = 0.1e1 / t33 ^ 2;
t60 = t58 * t48;
t42 = -t56 * t49 + t57 * t60;
t69 = t32 * t42;
t67 = t35 * t39;
t43 = t56 * t48 + t57 * t59;
t52 = 0.1e1 / t55 ^ 2;
t54 = 0.1e1 / t58 ^ 2;
t38 = 0.1e1 / (t43 ^ 2 * t54 * t52 + 0.1e1);
t51 = 0.1e1 / t55;
t66 = t38 * t51;
t65 = t42 ^ 2 * t32;
t46 = 0.1e1 / t48;
t64 = t46 * t51;
t63 = t52 * t57;
t53 = 0.1e1 / t58;
t47 = 0.1e1 / t48 ^ 2;
t41 = t49 * t61 - t60;
t36 = 0.1e1 / (t39 ^ 2 * t52 * t47 + 0.1e1);
t31 = 0.1e1 / t33;
t30 = (t39 * t46 * t63 + t56) * t36;
t29 = 0.1e1 / (0.1e1 + t65);
t28 = (t39 * t47 * t49 - t41 * t46) * t51 * t36;
t1 = [-t42 * t36 * t64, t30, 0, t28, 0, 0; (-t39 * t31 - (-t34 + (t64 * t67 + t34) * t36) * t65) * t29 (t30 * t67 * t69 + (-t58 * t55 * t31 - (t35 * t57 + (-t30 + t56) * t55 * t34) * t69) * t48) * t29, 0 (t43 * t31 - (t35 * t55 * t49 - t34 * t41 + (-t34 * t62 - t67) * t28) * t69) * t29, 0, 0; (t43 * t54 * t56 - t41 * t53) * t66 (-t43 * t53 * t63 - t49) * t38, 0, -t42 * t53 * t66, 0, 0;];
Ja_rot  = t1;
