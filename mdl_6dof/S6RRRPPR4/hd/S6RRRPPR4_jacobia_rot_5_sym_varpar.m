% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:08
% EndTime: 2019-02-26 22:05:08
% DurationCPUTime: 0.11s
% Computational Cost: add. (446->27), mult. (509->70), div. (100->11), fcn. (770->9), ass. (0->40)
t59 = sin(qJ(2));
t74 = t59 ^ 2;
t54 = qJ(3) + pkin(10);
t52 = sin(t54);
t53 = cos(t54);
t62 = cos(qJ(1));
t64 = t62 * t53;
t60 = sin(qJ(1));
t61 = cos(qJ(2));
t66 = t60 * t61;
t43 = t52 * t66 + t64;
t68 = t59 * t52;
t39 = atan2(-t43, t68);
t36 = sin(t39);
t37 = cos(t39);
t35 = -t36 * t43 + t37 * t68;
t34 = 0.1e1 / t35 ^ 2;
t65 = t62 * t52;
t46 = -t60 * t53 + t61 * t65;
t73 = t34 * t46;
t71 = t37 * t43;
t70 = t46 ^ 2 * t34;
t50 = 0.1e1 / t52;
t56 = 0.1e1 / t59;
t69 = t50 * t56;
t67 = t59 * t62;
t47 = t60 * t52 + t61 * t64;
t42 = 0.1e1 / t47 ^ 2;
t63 = t62 ^ 2 * t74 * t42;
t57 = 0.1e1 / t74;
t51 = 0.1e1 / t52 ^ 2;
t45 = t53 * t66 - t65;
t41 = 0.1e1 / t47;
t40 = 0.1e1 / (0.1e1 + t63);
t38 = 0.1e1 / (t43 ^ 2 * t57 * t51 + 0.1e1);
t33 = 0.1e1 / t35;
t32 = (t43 * t50 * t57 * t61 + t60) * t38;
t31 = 0.1e1 / (0.1e1 + t70);
t30 = (t43 * t51 * t53 - t45 * t50) * t56 * t38;
t1 = [-t46 * t38 * t69, t32, t30, 0, 0, 0; (-t43 * t33 - (-t36 + (t69 * t71 + t36) * t38) * t70) * t31 (t32 * t71 * t73 + (-t33 * t67 - (t37 * t61 + (-t32 + t60) * t59 * t36) * t73) * t52) * t31 (t47 * t33 - (t37 * t59 * t53 - t36 * t45 + (-t36 * t68 - t71) * t30) * t73) * t31, 0, 0, 0; (-t42 * t45 * t62 + t41 * t60) * t59 * t40 (-t41 * t61 * t62 - t53 * t63) * t40, -t46 * t42 * t40 * t67, 0, 0, 0;];
Ja_rot  = t1;
