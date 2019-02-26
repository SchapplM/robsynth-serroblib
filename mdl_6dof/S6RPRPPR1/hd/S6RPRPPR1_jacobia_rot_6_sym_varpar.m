% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:39:11
% EndTime: 2019-02-26 20:39:11
% DurationCPUTime: 0.11s
% Computational Cost: add. (395->23), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->36)
t55 = qJ(3) + pkin(10);
t52 = cos(t55);
t49 = sin(t55);
t56 = qJ(1) + pkin(9);
t50 = sin(t56);
t61 = t50 * t49;
t43 = atan2(-t61, -t52);
t41 = sin(t43);
t42 = cos(t43);
t34 = -t41 * t61 - t42 * t52;
t33 = 0.1e1 / t34 ^ 2;
t53 = cos(t56);
t66 = t33 * t53 ^ 2;
t54 = pkin(11) + qJ(6);
t48 = sin(t54);
t51 = cos(t54);
t58 = t53 * t51;
t40 = t50 * t48 + t52 * t58;
t38 = 0.1e1 / t40 ^ 2;
t59 = t53 * t48;
t39 = -t50 * t51 + t52 * t59;
t65 = t38 * t39;
t45 = t49 ^ 2;
t64 = t45 / t52 ^ 2;
t63 = t49 * t53;
t44 = 0.1e1 / (t50 ^ 2 * t64 + 0.1e1);
t62 = t50 * t44;
t60 = t50 * t52;
t57 = t39 ^ 2 * t38 + 0.1e1;
t46 = 0.1e1 / t52;
t37 = 0.1e1 / t40;
t36 = (0.1e1 + t64) * t62;
t35 = 0.1e1 / t57;
t32 = 0.1e1 / t34;
t31 = 0.1e1 / (t45 * t66 + 0.1e1);
t1 = [t46 * t44 * t63, 0, t36, 0, 0, 0; (-t32 * t61 - (-t42 * t45 * t46 * t62 + (t44 - 0.1e1) * t49 * t41) * t49 * t66) * t31, 0 (t52 * t32 - (-t41 * t60 + t42 * t49 + (t41 * t52 - t42 * t61) * t36) * t49 * t33) * t53 * t31, 0, 0, 0; ((-t48 * t60 - t58) * t37 - (-t51 * t60 + t59) * t65) * t35, 0 (-t37 * t48 + t51 * t65) * t35 * t63, 0, 0, t57 * t35;];
Ja_rot  = t1;
