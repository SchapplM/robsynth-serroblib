% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRP4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:31:55
% EndTime: 2019-02-26 20:31:55
% DurationCPUTime: 0.08s
% Computational Cost: add. (196->20), mult. (442->61), div. (52->9), fcn. (659->11), ass. (0->38)
t55 = sin(pkin(9));
t56 = cos(pkin(9));
t66 = sin(qJ(1));
t67 = cos(qJ(1));
t40 = -t66 * t55 - t67 * t56;
t68 = t40 ^ 2;
t53 = cos(qJ(4));
t41 = t67 * t55 - t66 * t56;
t51 = sin(qJ(4));
t60 = t41 * t51;
t38 = atan2(t60, t53);
t36 = sin(t38);
t37 = cos(t38);
t31 = t36 * t60 + t37 * t53;
t30 = 0.1e1 / t31 ^ 2;
t65 = t30 * t51;
t50 = sin(qJ(5));
t52 = cos(qJ(5));
t57 = t52 * t53;
t35 = -t40 * t57 + t41 * t50;
t33 = 0.1e1 / t35 ^ 2;
t58 = t50 * t53;
t34 = -t40 * t58 - t41 * t52;
t64 = t33 * t34;
t63 = t36 * t53;
t62 = t40 * t51;
t47 = t51 ^ 2;
t59 = t47 / t53 ^ 2;
t39 = 0.1e1 / (t41 ^ 2 * t59 + 0.1e1);
t61 = t41 * t39;
t54 = t34 ^ 2 * t33 + 0.1e1;
t48 = 0.1e1 / t53;
t32 = 0.1e1 / t35;
t29 = 0.1e1 / t31;
t28 = (0.1e1 + t59) * t61;
t27 = 0.1e1 / t54;
t26 = 0.1e1 / (t68 * t47 * t30 + 0.1e1);
t1 = [t48 * t39 * t62, 0, 0, t28, 0, 0; (t29 * t60 + (t37 * t47 * t48 * t61 + (-t39 + 0.1e1) * t51 * t36) * t68 * t65) * t26, 0, 0 (-t53 * t29 + (t41 * t63 - t37 * t51 + (t37 * t60 - t63) * t28) * t65) * t40 * t26, 0, 0; ((-t40 * t52 + t41 * t58) * t32 - (t40 * t50 + t41 * t57) * t64) * t27, 0, 0 (t32 * t50 - t52 * t64) * t27 * t62, t54 * t27, 0;];
Ja_rot  = t1;
