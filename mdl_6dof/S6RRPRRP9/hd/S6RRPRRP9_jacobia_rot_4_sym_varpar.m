% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP9_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:50:27
% EndTime: 2019-02-26 21:50:28
% DurationCPUTime: 0.09s
% Computational Cost: add. (221->25), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->41)
t58 = sin(qJ(2));
t60 = cos(qJ(2));
t61 = cos(qJ(1));
t59 = sin(qJ(1));
t65 = cos(pkin(6));
t63 = t59 * t65;
t49 = -t58 * t63 + t61 * t60;
t54 = pkin(11) + qJ(4);
t51 = sin(t54);
t52 = cos(t54);
t57 = sin(pkin(6));
t68 = t57 * t59;
t40 = t49 * t52 + t51 * t68;
t38 = 0.1e1 / t40 ^ 2;
t39 = t49 * t51 - t52 * t68;
t72 = t38 * t39;
t62 = t61 * t65;
t45 = t59 * t58 - t60 * t62;
t67 = t57 * t60;
t43 = atan2(-t45, -t67);
t42 = cos(t43);
t71 = t42 * t45;
t41 = sin(t43);
t35 = -t41 * t45 - t42 * t67;
t34 = 0.1e1 / t35 ^ 2;
t48 = t61 * t58 + t60 * t63;
t70 = t48 ^ 2 * t34;
t53 = 0.1e1 / t57;
t55 = 0.1e1 / t60;
t69 = t53 * t55;
t66 = t57 * t61;
t64 = t39 ^ 2 * t38 + 0.1e1;
t56 = 0.1e1 / t60 ^ 2;
t47 = t58 * t62 + t59 * t60;
t44 = 0.1e1 / (0.1e1 + t45 ^ 2 / t57 ^ 2 * t56);
t37 = 0.1e1 / t40;
t36 = 0.1e1 / t64;
t33 = 0.1e1 / t35;
t32 = 0.1e1 / (0.1e1 + t70);
t31 = (t45 * t56 * t58 + t47 * t55) * t53 * t44;
t1 = [t48 * t44 * t69, t31, 0, 0, 0, 0; (-t45 * t33 - (-t41 + (-t69 * t71 + t41) * t44) * t70) * t32 (t49 * t33 - (t42 * t57 * t58 - t41 * t47 + (t41 * t67 - t71) * t31) * t48 * t34) * t32, 0, 0, 0, 0; ((-t47 * t51 - t52 * t66) * t37 - (-t47 * t52 + t51 * t66) * t72) * t36 (-t51 * t37 + t52 * t72) * t48 * t36, 0, t64 * t36, 0, 0;];
Ja_rot  = t1;
