% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR7_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:32:06
% EndTime: 2019-02-26 21:32:06
% DurationCPUTime: 0.08s
% Computational Cost: add. (149->24), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->40)
t52 = sin(qJ(2));
t55 = cos(qJ(2));
t56 = cos(qJ(1));
t53 = sin(qJ(1));
t60 = cos(pkin(6));
t58 = t53 * t60;
t43 = t56 * t52 + t55 * t58;
t51 = sin(qJ(5));
t54 = cos(qJ(5));
t50 = sin(pkin(6));
t62 = t50 * t53;
t35 = t43 * t54 - t51 * t62;
t33 = 0.1e1 / t35 ^ 2;
t34 = t43 * t51 + t54 * t62;
t67 = t33 * t34;
t57 = t56 * t60;
t41 = t52 * t57 + t53 * t55;
t63 = t50 * t52;
t39 = atan2(-t41, t63);
t37 = cos(t39);
t66 = t37 * t41;
t36 = sin(t39);
t30 = -t36 * t41 + t37 * t63;
t29 = 0.1e1 / t30 ^ 2;
t44 = -t52 * t58 + t56 * t55;
t65 = t44 ^ 2 * t29;
t47 = 0.1e1 / t50;
t48 = 0.1e1 / t52;
t64 = t47 * t48;
t61 = t50 * t56;
t59 = t34 ^ 2 * t33 + 0.1e1;
t49 = 0.1e1 / t52 ^ 2;
t40 = t53 * t52 - t55 * t57;
t38 = 0.1e1 / (0.1e1 + t41 ^ 2 / t50 ^ 2 * t49);
t32 = 0.1e1 / t35;
t31 = 0.1e1 / t59;
t28 = 0.1e1 / t30;
t27 = 0.1e1 / (0.1e1 + t65);
t26 = (t41 * t49 * t55 + t40 * t48) * t47 * t38;
t1 = [-t44 * t38 * t64, t26, 0, 0, 0, 0; (-t41 * t28 - (-t36 + (t64 * t66 + t36) * t38) * t65) * t27 (-t43 * t28 - (t37 * t50 * t55 + t36 * t40 + (-t36 * t63 - t66) * t26) * t44 * t29) * t27, 0, 0, 0, 0; ((-t40 * t51 + t54 * t61) * t32 - (-t40 * t54 - t51 * t61) * t67) * t31 (t51 * t32 - t54 * t67) * t44 * t31, 0, 0, t59 * t31, 0;];
Ja_rot  = t1;
