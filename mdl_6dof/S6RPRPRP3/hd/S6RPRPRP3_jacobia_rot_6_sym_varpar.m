% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:44:51
% EndTime: 2019-02-26 20:44:51
% DurationCPUTime: 0.12s
% Computational Cost: add. (665->28), mult. (509->70), div. (100->11), fcn. (770->9), ass. (0->40)
t62 = sin(qJ(3));
t74 = t62 ^ 2;
t57 = pkin(10) + qJ(5);
t53 = sin(t57);
t55 = cos(t57);
t58 = qJ(1) + pkin(9);
t56 = cos(t58);
t54 = sin(t58);
t63 = cos(qJ(3));
t68 = t54 * t63;
t43 = t53 * t68 + t56 * t55;
t65 = t62 * t53;
t40 = atan2(-t43, t65);
t36 = sin(t40);
t37 = cos(t40);
t35 = -t36 * t43 + t37 * t65;
t34 = 0.1e1 / t35 ^ 2;
t66 = t56 * t63;
t46 = t53 * t66 - t54 * t55;
t73 = t34 * t46;
t71 = t37 * t43;
t70 = t46 ^ 2 * t34;
t50 = 0.1e1 / t53;
t60 = 0.1e1 / t62;
t69 = t50 * t60;
t67 = t56 * t62;
t47 = t54 * t53 + t55 * t66;
t42 = 0.1e1 / t47 ^ 2;
t64 = t56 ^ 2 * t74 * t42;
t61 = 0.1e1 / t74;
t51 = 0.1e1 / t53 ^ 2;
t45 = -t56 * t53 + t55 * t68;
t41 = 0.1e1 / t47;
t39 = 0.1e1 / (t43 ^ 2 * t61 * t51 + 0.1e1);
t38 = 0.1e1 / (0.1e1 + t64);
t33 = 0.1e1 / t35;
t32 = (t43 * t50 * t61 * t63 + t54) * t39;
t31 = 0.1e1 / (0.1e1 + t70);
t30 = (t43 * t51 * t55 - t45 * t50) * t60 * t39;
t1 = [-t46 * t39 * t69, 0, t32, 0, t30, 0; (-t43 * t33 - (-t36 + (t69 * t71 + t36) * t39) * t70) * t31, 0 (t32 * t71 * t73 + (-t33 * t67 - (t37 * t63 + (-t32 + t54) * t62 * t36) * t73) * t53) * t31, 0 (t47 * t33 - (t37 * t62 * t55 - t36 * t45 + (-t36 * t65 - t71) * t30) * t73) * t31, 0; (-t42 * t45 * t56 + t41 * t54) * t62 * t38, 0 (-t41 * t66 - t55 * t64) * t38, 0, -t46 * t42 * t38 * t67, 0;];
Ja_rot  = t1;
