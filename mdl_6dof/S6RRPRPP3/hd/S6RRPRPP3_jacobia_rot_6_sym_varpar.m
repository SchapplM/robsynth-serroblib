% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_rot = S6RRPRPP3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:56
% EndTime: 2019-02-26 21:35:56
% DurationCPUTime: 0.12s
% Computational Cost: add. (437->27), mult. (486->70), div. (108->11), fcn. (755->9), ass. (0->39)
t49 = pkin(9) + qJ(4);
t48 = cos(t49);
t47 = sin(t49);
t57 = cos(qJ(1));
t59 = t57 * t47;
t55 = sin(qJ(1));
t56 = cos(qJ(2));
t60 = t55 * t56;
t40 = t48 * t60 - t59;
t54 = sin(qJ(2));
t61 = t54 * t48;
t37 = atan2(-t40, t61);
t34 = sin(t37);
t35 = cos(t37);
t33 = -t34 * t40 + t35 * t61;
t32 = 0.1e1 / t33 ^ 2;
t58 = t57 * t48;
t43 = t55 * t47 + t56 * t58;
t68 = t32 * t43;
t66 = t35 * t40;
t42 = t55 * t48 - t56 * t59;
t51 = 0.1e1 / t54 ^ 2;
t53 = 0.1e1 / t57 ^ 2;
t38 = 0.1e1 / (t42 ^ 2 * t53 * t51 + 0.1e1);
t50 = 0.1e1 / t54;
t65 = t38 * t50;
t64 = t43 ^ 2 * t32;
t45 = 0.1e1 / t48;
t63 = t45 * t50;
t62 = t51 * t56;
t52 = 0.1e1 / t57;
t46 = 0.1e1 / t48 ^ 2;
t39 = t47 * t60 + t58;
t36 = 0.1e1 / (t40 ^ 2 * t51 * t46 + 0.1e1);
t31 = 0.1e1 / t33;
t30 = (t40 * t45 * t62 + t55) * t36;
t29 = 0.1e1 / (0.1e1 + t64);
t28 = (-t40 * t46 * t47 + t39 * t45) * t50 * t36;
t1 = [-t43 * t36 * t63, t30, 0, t28, 0, 0; (-t40 * t31 - (-t34 + (t63 * t66 + t34) * t36) * t64) * t29 (t30 * t66 * t68 + (-t57 * t54 * t31 - (t35 * t56 + (-t30 + t55) * t54 * t34) * t68) * t48) * t29, 0 (t42 * t31 - (-t35 * t54 * t47 + t34 * t39 + (-t34 * t61 - t66) * t28) * t68) * t29, 0, 0; (t42 * t53 * t55 + t39 * t52) * t65 (-t42 * t52 * t62 + t47) * t38, 0, -t43 * t52 * t65, 0, 0;];
Ja_rot  = t1;
