% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR10_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:03
% EndTime: 2019-02-26 21:06:03
% DurationCPUTime: 0.10s
% Computational Cost: add. (160->24), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->40)
t54 = cos(qJ(3));
t68 = t54 ^ 2;
t51 = sin(qJ(3));
t50 = sin(qJ(4));
t55 = cos(qJ(1));
t58 = t55 * t50;
t52 = sin(qJ(1));
t53 = cos(qJ(4));
t61 = t52 * t53;
t41 = t51 * t58 + t61;
t59 = t54 * t50;
t36 = atan2(t41, t59);
t32 = sin(t36);
t33 = cos(t36);
t31 = t32 * t41 + t33 * t59;
t30 = 0.1e1 / t31 ^ 2;
t57 = t55 * t53;
t62 = t52 * t50;
t39 = t51 * t62 - t57;
t67 = t30 * t39;
t65 = t33 * t41;
t64 = t39 ^ 2 * t30;
t44 = 0.1e1 / t50;
t48 = 0.1e1 / t54;
t63 = t44 * t48;
t60 = t52 * t54;
t40 = t51 * t61 + t58;
t38 = 0.1e1 / t40 ^ 2;
t56 = t52 ^ 2 * t68 * t38;
t49 = 0.1e1 / t68;
t45 = 0.1e1 / t50 ^ 2;
t42 = t51 * t57 - t62;
t37 = 0.1e1 / t40;
t35 = 0.1e1 / (t41 ^ 2 * t49 * t45 + 0.1e1);
t34 = 0.1e1 / (0.1e1 + t56);
t29 = 0.1e1 / t31;
t28 = (t41 * t44 * t49 * t51 + t55) * t35;
t27 = 0.1e1 / (0.1e1 + t64);
t26 = (-t41 * t45 * t53 + t42 * t44) * t48 * t35;
t1 = [-t39 * t35 * t63, 0, t28, t26, 0, 0; (t41 * t29 - (-t32 + (-t63 * t65 + t32) * t35) * t64) * t27, 0 (-t28 * t65 * t67 + (t29 * t60 - (-t33 * t51 + (-t28 + t55) * t54 * t32) * t67) * t50) * t27 (t40 * t29 - (t33 * t54 * t53 + t32 * t42 + (-t32 * t59 + t65) * t26) * t67) * t27, 0, 0; (-t38 * t42 * t52 + t37 * t55) * t54 * t34, 0 (-t37 * t51 * t52 - t53 * t56) * t34, t39 * t38 * t34 * t60, 0, 0;];
Ja_rot  = t1;
