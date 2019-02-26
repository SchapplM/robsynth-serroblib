% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRRRR7
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR7_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_jacobia_rot_3_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:41
% EndTime: 2019-02-26 22:50:41
% DurationCPUTime: 0.11s
% Computational Cost: add. (173->24), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->40)
t53 = sin(qJ(2));
t56 = cos(qJ(2));
t57 = cos(qJ(1));
t54 = sin(qJ(1));
t61 = cos(pkin(6));
t59 = t54 * t61;
t46 = -t53 * t59 + t57 * t56;
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t51 = sin(pkin(6));
t64 = t51 * t54;
t37 = t46 * t55 + t52 * t64;
t35 = 0.1e1 / t37 ^ 2;
t36 = t46 * t52 - t55 * t64;
t68 = t35 * t36;
t58 = t57 * t61;
t42 = t54 * t53 - t56 * t58;
t63 = t51 * t56;
t40 = atan2(-t42, -t63);
t39 = cos(t40);
t67 = t39 * t42;
t38 = sin(t40);
t32 = -t38 * t42 - t39 * t63;
t31 = 0.1e1 / t32 ^ 2;
t45 = t57 * t53 + t56 * t59;
t66 = t45 ^ 2 * t31;
t48 = 0.1e1 / t51;
t49 = 0.1e1 / t56;
t65 = t48 * t49;
t62 = t51 * t57;
t60 = t36 ^ 2 * t35 + 0.1e1;
t50 = 0.1e1 / t56 ^ 2;
t44 = t53 * t58 + t54 * t56;
t41 = 0.1e1 / (0.1e1 + t42 ^ 2 / t51 ^ 2 * t50);
t34 = 0.1e1 / t37;
t33 = 0.1e1 / t60;
t30 = 0.1e1 / t32;
t29 = 0.1e1 / (0.1e1 + t66);
t28 = (t42 * t50 * t53 + t44 * t49) * t48 * t41;
t1 = [t45 * t41 * t65, t28, 0, 0, 0, 0; (-t42 * t30 - (-t38 + (-t65 * t67 + t38) * t41) * t66) * t29 (t46 * t30 - (t39 * t51 * t53 - t38 * t44 + (t38 * t63 - t67) * t28) * t45 * t31) * t29, 0, 0, 0, 0; ((-t44 * t52 - t55 * t62) * t34 - (-t44 * t55 + t52 * t62) * t68) * t33 (-t52 * t34 + t55 * t68) * t45 * t33, t60 * t33, 0, 0, 0;];
Ja_rot  = t1;
