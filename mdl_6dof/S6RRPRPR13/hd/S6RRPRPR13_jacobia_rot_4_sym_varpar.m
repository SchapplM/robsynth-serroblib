% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR13_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:45:00
% EndTime: 2019-02-26 21:45:00
% DurationCPUTime: 0.08s
% Computational Cost: add. (149->24), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->40)
t54 = sin(qJ(2));
t57 = cos(qJ(2));
t58 = cos(qJ(1));
t55 = sin(qJ(1));
t62 = cos(pkin(6));
t60 = t55 * t62;
t45 = t58 * t54 + t57 * t60;
t53 = sin(qJ(4));
t56 = cos(qJ(4));
t52 = sin(pkin(6));
t64 = t52 * t55;
t37 = t45 * t53 + t56 * t64;
t35 = 0.1e1 / t37 ^ 2;
t36 = -t45 * t56 + t53 * t64;
t69 = t35 * t36;
t59 = t58 * t62;
t43 = t54 * t59 + t55 * t57;
t65 = t52 * t54;
t41 = atan2(-t43, t65);
t39 = cos(t41);
t68 = t39 * t43;
t38 = sin(t41);
t32 = -t38 * t43 + t39 * t65;
t31 = 0.1e1 / t32 ^ 2;
t46 = -t54 * t60 + t58 * t57;
t67 = t46 ^ 2 * t31;
t49 = 0.1e1 / t52;
t50 = 0.1e1 / t54;
t66 = t49 * t50;
t63 = t52 * t58;
t61 = t36 ^ 2 * t35 + 0.1e1;
t51 = 0.1e1 / t54 ^ 2;
t42 = t55 * t54 - t57 * t59;
t40 = 0.1e1 / (0.1e1 + t43 ^ 2 / t52 ^ 2 * t51);
t34 = 0.1e1 / t37;
t33 = 0.1e1 / t61;
t30 = 0.1e1 / t32;
t29 = 0.1e1 / (0.1e1 + t67);
t28 = (t43 * t51 * t57 + t42 * t50) * t49 * t40;
t1 = [-t46 * t40 * t66, t28, 0, 0, 0, 0; (-t43 * t30 - (-t38 + (t66 * t68 + t38) * t40) * t67) * t29 (-t45 * t30 - (t39 * t52 * t57 + t38 * t42 + (-t38 * t65 - t68) * t28) * t46 * t31) * t29, 0, 0, 0, 0; ((t42 * t56 + t53 * t63) * t34 - (-t42 * t53 + t56 * t63) * t69) * t33 (-t56 * t34 - t53 * t69) * t46 * t33, 0, t61 * t33, 0, 0;];
Ja_rot  = t1;
