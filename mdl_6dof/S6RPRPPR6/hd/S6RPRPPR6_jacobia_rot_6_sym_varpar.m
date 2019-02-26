% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:41:59
% EndTime: 2019-02-26 20:41:59
% DurationCPUTime: 0.10s
% Computational Cost: add. (263->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
t53 = qJ(3) + pkin(9);
t49 = sin(t53);
t51 = cos(t53);
t55 = cos(qJ(1));
t57 = t55 * t51;
t43 = atan2(-t57, t49);
t41 = sin(t43);
t42 = cos(t43);
t34 = -t41 * t57 + t42 * t49;
t33 = 0.1e1 / t34 ^ 2;
t54 = sin(qJ(1));
t67 = t33 * t54 ^ 2;
t52 = pkin(10) + qJ(6);
t48 = sin(t52);
t59 = t55 * t48;
t50 = cos(t52);
t61 = t54 * t50;
t40 = t49 * t61 + t59;
t38 = 0.1e1 / t40 ^ 2;
t58 = t55 * t50;
t62 = t54 * t48;
t39 = t49 * t62 - t58;
t66 = t38 * t39;
t65 = t41 * t49;
t47 = t51 ^ 2;
t64 = 0.1e1 / t49 ^ 2 * t47;
t63 = t51 * t54;
t44 = 0.1e1 / (t55 ^ 2 * t64 + 0.1e1);
t60 = t55 * t44;
t56 = t39 ^ 2 * t38 + 0.1e1;
t45 = 0.1e1 / t49;
t37 = 0.1e1 / t40;
t36 = (0.1e1 + t64) * t60;
t35 = 0.1e1 / t56;
t32 = 0.1e1 / t34;
t31 = 0.1e1 / (t47 * t67 + 0.1e1);
t1 = [t45 * t44 * t63, 0, t36, 0, 0, 0; (-t32 * t57 + (-t42 * t45 * t47 * t60 + (-t44 + 0.1e1) * t51 * t41) * t51 * t67) * t31, 0 (t49 * t32 + (t55 * t65 + t42 * t51 + (-t42 * t57 - t65) * t36) * t51 * t33) * t54 * t31, 0, 0, 0; ((t49 * t59 + t61) * t37 - (t49 * t58 - t62) * t66) * t35, 0 (t37 * t48 - t50 * t66) * t35 * t63, 0, 0, t56 * t35;];
Ja_rot  = t1;
