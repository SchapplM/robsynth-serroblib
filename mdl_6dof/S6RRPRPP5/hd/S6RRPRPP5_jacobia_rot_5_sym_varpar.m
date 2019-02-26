% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPP5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_jacobia_rot_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:06
% EndTime: 2019-02-26 21:37:06
% DurationCPUTime: 0.10s
% Computational Cost: add. (159->25), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->40)
t53 = cos(qJ(2));
t67 = t53 ^ 2;
t50 = sin(qJ(2));
t49 = sin(qJ(4));
t54 = cos(qJ(1));
t57 = t54 * t49;
t51 = sin(qJ(1));
t52 = cos(qJ(4));
t60 = t51 * t52;
t40 = t50 * t60 + t57;
t59 = t53 * t52;
t35 = atan2(t40, t59);
t31 = sin(t35);
t32 = cos(t35);
t30 = t31 * t40 + t32 * t59;
t29 = 0.1e1 / t30 ^ 2;
t56 = t54 * t52;
t61 = t51 * t49;
t38 = -t50 * t56 + t61;
t66 = t29 * t38;
t64 = t32 * t40;
t63 = t38 ^ 2 * t29;
t43 = 0.1e1 / t52;
t46 = 0.1e1 / t53;
t62 = t43 * t46;
t58 = t53 * t54;
t39 = t50 * t57 + t60;
t37 = 0.1e1 / t39 ^ 2;
t55 = t54 ^ 2 * t67 * t37;
t47 = 0.1e1 / t67;
t44 = 0.1e1 / t52 ^ 2;
t41 = -t50 * t61 + t56;
t36 = 0.1e1 / t39;
t34 = 0.1e1 / (t40 ^ 2 * t47 * t44 + 0.1e1);
t33 = 0.1e1 / (0.1e1 + t55);
t28 = 0.1e1 / t30;
t27 = (t40 * t43 * t47 * t50 + t51) * t34;
t26 = 0.1e1 / (0.1e1 + t63);
t25 = (t40 * t44 * t49 + t41 * t43) * t46 * t34;
t1 = [-t38 * t34 * t62, t27, 0, t25, 0, 0; (t40 * t28 - (-t31 + (-t62 * t64 + t31) * t34) * t63) * t26 (-t27 * t64 * t66 + (-t28 * t58 - (-t32 * t50 + (-t27 + t51) * t53 * t31) * t66) * t52) * t26, 0 (t39 * t28 - (-t32 * t53 * t49 + t31 * t41 + (-t31 * t59 + t64) * t25) * t66) * t26, 0, 0; (t37 * t41 * t54 + t36 * t51) * t53 * t33 (t36 * t50 * t54 + t49 * t55) * t33, 0, -t38 * t37 * t33 * t58, 0, 0;];
Ja_rot  = t1;
