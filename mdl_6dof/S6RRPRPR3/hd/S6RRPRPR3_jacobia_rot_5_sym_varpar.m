% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:38:51
% EndTime: 2019-02-26 21:38:51
% DurationCPUTime: 0.10s
% Computational Cost: add. (286->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
t54 = qJ(2) + pkin(10);
t52 = cos(t54);
t50 = sin(t54);
t55 = sin(qJ(1));
t61 = t55 * t50;
t44 = atan2(-t61, -t52);
t42 = sin(t44);
t43 = cos(t44);
t35 = -t42 * t61 - t43 * t52;
t34 = 0.1e1 / t35 ^ 2;
t56 = cos(qJ(1));
t68 = t34 * t56 ^ 2;
t53 = qJ(4) + pkin(11);
t51 = cos(t53);
t58 = t56 * t51;
t49 = sin(t53);
t62 = t55 * t49;
t41 = t52 * t58 + t62;
t39 = 0.1e1 / t41 ^ 2;
t59 = t56 * t49;
t60 = t55 * t51;
t40 = t52 * t59 - t60;
t67 = t39 * t40;
t66 = t42 * t52;
t46 = t50 ^ 2;
t65 = t46 / t52 ^ 2;
t64 = t50 * t56;
t45 = 0.1e1 / (t55 ^ 2 * t65 + 0.1e1);
t63 = t55 * t45;
t57 = t40 ^ 2 * t39 + 0.1e1;
t47 = 0.1e1 / t52;
t38 = 0.1e1 / t41;
t37 = (0.1e1 + t65) * t63;
t36 = 0.1e1 / t57;
t33 = 0.1e1 / t35;
t32 = 0.1e1 / (t46 * t68 + 0.1e1);
t1 = [t47 * t45 * t64, t37, 0, 0, 0, 0; (-t33 * t61 - (-t43 * t46 * t47 * t63 + (t45 - 0.1e1) * t50 * t42) * t50 * t68) * t32 (t52 * t33 - (-t55 * t66 + t43 * t50 + (-t43 * t61 + t66) * t37) * t50 * t34) * t56 * t32, 0, 0, 0, 0; ((-t52 * t62 - t58) * t38 - (-t52 * t60 + t59) * t67) * t36 (-t38 * t49 + t51 * t67) * t36 * t64, 0, t57 * t36, 0, 0;];
Ja_rot  = t1;
