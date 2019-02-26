% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6RRRPPR7_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:11
% EndTime: 2019-02-26 22:07:11
% DurationCPUTime: 0.11s
% Computational Cost: add. (117->25), mult. (363->67), div. (53->9), fcn. (527->11), ass. (0->39)
t57 = sin(qJ(1));
t58 = cos(qJ(3));
t59 = cos(qJ(2));
t55 = sin(qJ(3));
t60 = cos(qJ(1));
t62 = t60 * t55;
t43 = -t57 * t58 + t59 * t62;
t61 = t60 * t58;
t44 = t57 * t55 + t59 * t61;
t53 = sin(pkin(10));
t54 = cos(pkin(10));
t35 = t43 * t53 + t44 * t54;
t33 = 0.1e1 / t35 ^ 2;
t34 = -t43 * t54 + t44 * t53;
t70 = t33 * t34;
t69 = t34 ^ 2 * t33;
t56 = sin(qJ(2));
t64 = t57 * t56;
t48 = atan2(t64, t59);
t45 = sin(t48);
t46 = cos(t48);
t38 = t45 * t64 + t46 * t59;
t37 = 0.1e1 / t38 ^ 2;
t68 = t37 * t60 ^ 2;
t50 = t56 ^ 2;
t67 = t50 / t59 ^ 2;
t66 = t56 * t60;
t47 = 0.1e1 / (t57 ^ 2 * t67 + 0.1e1);
t65 = t57 * t47;
t63 = t57 * t59;
t51 = 0.1e1 / t59;
t42 = -t58 * t63 + t62;
t41 = -t55 * t63 - t61;
t39 = (0.1e1 + t67) * t65;
t36 = 0.1e1 / t38;
t32 = 0.1e1 / t35;
t31 = 0.1e1 / (t50 * t68 + 0.1e1);
t30 = 0.1e1 / (0.1e1 + t69);
t1 = [t51 * t47 * t66, t39, 0, 0, 0, 0; (t36 * t64 + (t46 * t50 * t51 * t65 + (-t47 + 0.1e1) * t56 * t45) * t56 * t68) * t31 (-t59 * t36 + (t45 * t63 - t46 * t56 + (-t45 * t59 + t46 * t64) * t39) * t56 * t37) * t60 * t31, 0, 0, 0, 0; ((-t41 * t54 + t42 * t53) * t32 - (t41 * t53 + t42 * t54) * t70) * t30 ((-t53 * t58 + t54 * t55) * t32 - (-t53 * t55 - t54 * t58) * t70) * t30 * t66 (-t35 * t32 - t69) * t30, 0, 0, 0;];
Ja_rot  = t1;
