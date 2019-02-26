% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR7_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:37
% EndTime: 2019-02-26 20:07:37
% DurationCPUTime: 0.09s
% Computational Cost: add. (253->26), mult. (731->68), div. (57->9), fcn. (1053->11), ass. (0->42)
t53 = sin(pkin(11));
t55 = cos(pkin(11));
t60 = cos(qJ(2));
t56 = cos(pkin(6));
t58 = sin(qJ(2));
t63 = t56 * t58;
t48 = t53 * t60 + t55 * t63;
t57 = sin(qJ(3));
t54 = sin(pkin(6));
t59 = cos(qJ(3));
t65 = t54 * t59;
t38 = t48 * t57 + t55 * t65;
t66 = t54 * t57;
t51 = -t56 * t59 + t58 * t66;
t37 = atan2(-t38, t51);
t33 = sin(t37);
t34 = cos(t37);
t32 = -t33 * t38 + t34 * t51;
t31 = 0.1e1 / t32 ^ 2;
t50 = -t53 * t63 + t55 * t60;
t41 = t50 * t57 - t53 * t65;
t68 = t31 * t41;
t46 = 0.1e1 / t51 ^ 2;
t67 = t38 * t46;
t64 = t54 * t60;
t62 = t56 * t60;
t61 = -t33 * t51 - t34 * t38;
t52 = t56 * t57 + t58 * t65;
t49 = t53 * t62 + t55 * t58;
t47 = -t53 * t58 + t55 * t62;
t45 = 0.1e1 / t51;
t44 = 0.1e1 / t49 ^ 2;
t43 = 0.1e1 / t49;
t42 = t50 * t59 + t53 * t66;
t40 = t48 * t59 - t55 * t66;
t36 = 0.1e1 / (t38 ^ 2 * t46 + 0.1e1);
t35 = 0.1e1 / (t42 ^ 2 * t44 + 0.1e1);
t30 = 0.1e1 / t32;
t29 = 0.1e1 / (t41 ^ 2 * t31 + 0.1e1);
t28 = (-t45 * t47 + t64 * t67) * t57 * t36;
t27 = (-t40 * t45 + t52 * t67) * t36;
t1 = [0, t28, t27, 0, 0, 0; 0 (-t49 * t57 * t30 - ((-t33 * t47 + t34 * t64) * t57 + t61 * t28) * t68) * t29 (t42 * t30 - (t61 * t27 - t33 * t40 + t34 * t52) * t68) * t29, 0, 0, 0; 0 (-t42 * t44 * t50 - t43 * t49 * t59) * t35, -t41 * t43 * t35, 0, 0, 0;];
Ja_rot  = t1;
