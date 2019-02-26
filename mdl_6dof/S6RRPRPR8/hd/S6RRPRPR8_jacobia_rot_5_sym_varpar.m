% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR8_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:47
% EndTime: 2019-02-26 21:41:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (446->27), mult. (509->70), div. (100->11), fcn. (770->9), ass. (0->40)
t57 = sin(qJ(2));
t72 = t57 ^ 2;
t52 = pkin(10) + qJ(4);
t50 = sin(t52);
t51 = cos(t52);
t60 = cos(qJ(1));
t62 = t60 * t51;
t58 = sin(qJ(1));
t59 = cos(qJ(2));
t64 = t58 * t59;
t41 = t50 * t64 + t62;
t66 = t57 * t50;
t37 = atan2(-t41, t66);
t34 = sin(t37);
t35 = cos(t37);
t33 = -t34 * t41 + t35 * t66;
t32 = 0.1e1 / t33 ^ 2;
t63 = t60 * t50;
t44 = -t58 * t51 + t59 * t63;
t71 = t32 * t44;
t69 = t35 * t41;
t68 = t44 ^ 2 * t32;
t48 = 0.1e1 / t50;
t54 = 0.1e1 / t57;
t67 = t48 * t54;
t65 = t57 * t60;
t45 = t58 * t50 + t59 * t62;
t40 = 0.1e1 / t45 ^ 2;
t61 = t60 ^ 2 * t72 * t40;
t55 = 0.1e1 / t72;
t49 = 0.1e1 / t50 ^ 2;
t43 = t51 * t64 - t63;
t39 = 0.1e1 / t45;
t38 = 0.1e1 / (0.1e1 + t61);
t36 = 0.1e1 / (t41 ^ 2 * t55 * t49 + 0.1e1);
t31 = 0.1e1 / t33;
t30 = (t41 * t48 * t55 * t59 + t58) * t36;
t29 = 0.1e1 / (0.1e1 + t68);
t28 = (t41 * t49 * t51 - t43 * t48) * t54 * t36;
t1 = [-t44 * t36 * t67, t30, 0, t28, 0, 0; (-t41 * t31 - (-t34 + (t67 * t69 + t34) * t36) * t68) * t29 (t30 * t69 * t71 + (-t31 * t65 - (t35 * t59 + (-t30 + t58) * t57 * t34) * t71) * t50) * t29, 0 (t45 * t31 - (t35 * t57 * t51 - t34 * t43 + (-t34 * t66 - t69) * t28) * t71) * t29, 0, 0; (-t40 * t43 * t60 + t39 * t58) * t57 * t38 (-t39 * t59 * t60 - t51 * t61) * t38, 0, -t44 * t40 * t38 * t65, 0, 0;];
Ja_rot  = t1;
