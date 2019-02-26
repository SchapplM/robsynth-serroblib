% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPP2
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
% Datum: 2019-02-26 21:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPP2_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_jacobia_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:21
% EndTime: 2019-02-26 21:35:21
% DurationCPUTime: 0.07s
% Computational Cost: add. (238->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
t46 = qJ(2) + pkin(9);
t45 = cos(t46);
t44 = sin(t46);
t48 = sin(qJ(1));
t56 = t48 * t44;
t39 = atan2(-t56, -t45);
t35 = sin(t39);
t36 = cos(t39);
t30 = -t35 * t56 - t36 * t45;
t29 = 0.1e1 / t30 ^ 2;
t50 = cos(qJ(1));
t62 = t29 * t50 ^ 2;
t49 = cos(qJ(4));
t52 = t50 * t49;
t47 = sin(qJ(4));
t55 = t48 * t47;
t38 = t45 * t52 + t55;
t34 = 0.1e1 / t38 ^ 2;
t53 = t50 * t47;
t54 = t48 * t49;
t37 = t45 * t53 - t54;
t61 = t34 * t37;
t60 = t35 * t45;
t41 = t44 ^ 2;
t59 = t41 / t45 ^ 2;
t58 = t44 * t50;
t40 = 0.1e1 / (t48 ^ 2 * t59 + 0.1e1);
t57 = t48 * t40;
t51 = t37 ^ 2 * t34 + 0.1e1;
t42 = 0.1e1 / t45;
t33 = 0.1e1 / t38;
t32 = 0.1e1 / t51;
t31 = (0.1e1 + t59) * t57;
t28 = 0.1e1 / t30;
t27 = 0.1e1 / (t41 * t62 + 0.1e1);
t1 = [t42 * t40 * t58, t31, 0, 0, 0, 0; (-t28 * t56 - (-t36 * t41 * t42 * t57 + (t40 - 0.1e1) * t44 * t35) * t44 * t62) * t27 (t45 * t28 - (-t48 * t60 + t36 * t44 + (-t36 * t56 + t60) * t31) * t44 * t29) * t50 * t27, 0, 0, 0, 0; ((-t45 * t55 - t52) * t33 - (-t45 * t54 + t53) * t61) * t32 (-t33 * t47 + t49 * t61) * t32 * t58, 0, t51 * t32, 0, 0;];
Ja_rot  = t1;
