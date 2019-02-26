% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRP2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:30:41
% EndTime: 2019-02-26 20:30:41
% DurationCPUTime: 0.11s
% Computational Cost: add. (572->28), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->42)
t57 = pkin(10) + qJ(4);
t53 = sin(t57);
t75 = t53 ^ 2;
t55 = cos(t57);
t58 = qJ(1) + pkin(9);
t56 = cos(t58);
t62 = cos(qJ(5));
t64 = t56 * t62;
t54 = sin(t58);
t61 = sin(qJ(5));
t67 = t54 * t61;
t42 = t55 * t67 + t64;
t68 = t53 * t61;
t39 = atan2(-t42, t68);
t36 = sin(t39);
t37 = cos(t39);
t34 = -t36 * t42 + t37 * t68;
t33 = 0.1e1 / t34 ^ 2;
t65 = t56 * t61;
t66 = t54 * t62;
t45 = t55 * t65 - t66;
t74 = t33 * t45;
t72 = t37 * t42;
t71 = t45 ^ 2 * t33;
t50 = 0.1e1 / t53;
t59 = 0.1e1 / t61;
t70 = t50 * t59;
t69 = t53 * t56;
t46 = t55 * t64 + t67;
t41 = 0.1e1 / t46 ^ 2;
t63 = t56 ^ 2 * t75 * t41;
t60 = 0.1e1 / t61 ^ 2;
t51 = 0.1e1 / t75;
t44 = t55 * t66 - t65;
t40 = 0.1e1 / t46;
t38 = 0.1e1 / (t42 ^ 2 * t51 * t60 + 0.1e1);
t35 = 0.1e1 / (0.1e1 + t63);
t32 = 0.1e1 / t34;
t31 = (t42 * t51 * t55 * t59 + t54) * t38;
t30 = 0.1e1 / (0.1e1 + t71);
t29 = (t42 * t60 * t62 - t44 * t59) * t50 * t38;
t1 = [-t45 * t38 * t70, 0, 0, t31, t29, 0; (-t42 * t32 - (-t36 + (t70 * t72 + t36) * t38) * t71) * t30, 0, 0 (t31 * t72 * t74 + (-t32 * t69 - (t37 * t55 + (-t31 + t54) * t53 * t36) * t74) * t61) * t30 (t46 * t32 - (t37 * t53 * t62 - t36 * t44 + (-t36 * t68 - t72) * t29) * t74) * t30, 0; (-t41 * t44 * t56 + t40 * t54) * t53 * t35, 0, 0 (-t40 * t55 * t56 - t62 * t63) * t35, -t45 * t41 * t35 * t69, 0;];
Ja_rot  = t1;
