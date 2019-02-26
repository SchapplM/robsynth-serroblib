% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:35:57
% EndTime: 2019-02-26 20:35:57
% DurationCPUTime: 0.11s
% Computational Cost: add. (281->21), mult. (251->57), div. (57->9), fcn. (367->9), ass. (0->37)
t53 = qJ(1) + pkin(10);
t49 = sin(t53);
t70 = t49 ^ 2;
t58 = sin(qJ(4));
t50 = cos(t53);
t59 = cos(qJ(4));
t64 = t50 * t59;
t47 = atan2(-t64, t58);
t45 = sin(t47);
t46 = cos(t47);
t39 = -t45 * t64 + t46 * t58;
t38 = 0.1e1 / t39 ^ 2;
t69 = t38 * t59;
t57 = qJ(5) + qJ(6);
t51 = sin(t57);
t52 = cos(t57);
t62 = t52 * t58;
t44 = t49 * t62 + t50 * t51;
t42 = 0.1e1 / t44 ^ 2;
t63 = t51 * t58;
t43 = t49 * t63 - t50 * t52;
t68 = t42 * t43;
t67 = t45 * t58;
t66 = t49 * t59;
t56 = t59 ^ 2;
t61 = 0.1e1 / t58 ^ 2 * t56;
t48 = 0.1e1 / (t50 ^ 2 * t61 + 0.1e1);
t65 = t50 * t48;
t60 = t43 ^ 2 * t42 + 0.1e1;
t54 = 0.1e1 / t58;
t41 = 0.1e1 / t44;
t40 = (0.1e1 + t61) * t65;
t37 = 0.1e1 / t39;
t36 = 0.1e1 / t60;
t35 = 0.1e1 / (t70 * t56 * t38 + 0.1e1);
t34 = t60 * t36;
t1 = [t54 * t48 * t66, 0, 0, t40, 0, 0; (-t37 * t64 + (-t46 * t54 * t56 * t65 + (-t48 + 0.1e1) * t59 * t45) * t70 * t69) * t35, 0, 0 (t58 * t37 + (t50 * t67 + t46 * t59 + (-t46 * t64 - t67) * t40) * t69) * t49 * t35, 0, 0; ((t49 * t52 + t50 * t63) * t41 - (-t49 * t51 + t50 * t62) * t68) * t36, 0, 0 (t41 * t51 - t52 * t68) * t36 * t66, t34, t34;];
Ja_rot  = t1;
