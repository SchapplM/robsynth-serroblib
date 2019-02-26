% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR7_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobia_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:01
% EndTime: 2019-02-26 22:07:01
% DurationCPUTime: 0.10s
% Computational Cost: add. (159->26), mult. (509->70), div. (100->11), fcn. (770->9), ass. (0->39)
t51 = sin(qJ(2));
t67 = t51 ^ 2;
t50 = sin(qJ(3));
t53 = cos(qJ(3));
t55 = cos(qJ(1));
t57 = t55 * t53;
t52 = sin(qJ(1));
t54 = cos(qJ(2));
t59 = t52 * t54;
t37 = t50 * t59 + t57;
t61 = t51 * t50;
t34 = atan2(-t37, t61);
t30 = sin(t34);
t31 = cos(t34);
t29 = -t30 * t37 + t31 * t61;
t28 = 0.1e1 / t29 ^ 2;
t58 = t55 * t50;
t40 = -t52 * t53 + t54 * t58;
t66 = t28 * t40;
t64 = t31 * t37;
t63 = t40 ^ 2 * t28;
t44 = 0.1e1 / t50;
t47 = 0.1e1 / t51;
t62 = t44 * t47;
t60 = t51 * t55;
t41 = t52 * t50 + t54 * t57;
t36 = 0.1e1 / t41 ^ 2;
t56 = t55 ^ 2 * t67 * t36;
t48 = 0.1e1 / t67;
t45 = 0.1e1 / t50 ^ 2;
t39 = t53 * t59 - t58;
t35 = 0.1e1 / t41;
t33 = 0.1e1 / (t37 ^ 2 * t48 * t45 + 0.1e1);
t32 = 0.1e1 / (0.1e1 + t56);
t27 = 0.1e1 / t29;
t26 = (t37 * t44 * t48 * t54 + t52) * t33;
t25 = 0.1e1 / (0.1e1 + t63);
t24 = (t37 * t45 * t53 - t39 * t44) * t47 * t33;
t1 = [-t40 * t33 * t62, t26, t24, 0, 0, 0; (-t37 * t27 - (-t30 + (t62 * t64 + t30) * t33) * t63) * t25 (t26 * t64 * t66 + (-t27 * t60 - (t31 * t54 + (-t26 + t52) * t51 * t30) * t66) * t50) * t25 (t41 * t27 - (t31 * t51 * t53 - t30 * t39 + (-t30 * t61 - t64) * t24) * t66) * t25, 0, 0, 0; (-t36 * t39 * t55 + t35 * t52) * t51 * t32 (-t35 * t54 * t55 - t53 * t56) * t32, -t40 * t36 * t32 * t60, 0, 0, 0;];
Ja_rot  = t1;
