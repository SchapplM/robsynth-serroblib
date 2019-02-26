% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP1
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

function Ja_rot = S6RPPRRP1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:30:10
% EndTime: 2019-02-26 20:30:10
% DurationCPUTime: 0.10s
% Computational Cost: add. (347->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
t50 = pkin(10) + qJ(4);
t48 = cos(t50);
t46 = sin(t50);
t51 = qJ(1) + pkin(9);
t47 = sin(t51);
t59 = t47 * t46;
t41 = atan2(-t59, -t48);
t39 = sin(t41);
t40 = cos(t41);
t32 = -t39 * t59 - t40 * t48;
t31 = 0.1e1 / t32 ^ 2;
t49 = cos(t51);
t65 = t31 * t49 ^ 2;
t53 = cos(qJ(5));
t55 = t49 * t53;
t52 = sin(qJ(5));
t58 = t47 * t52;
t38 = t48 * t55 + t58;
t36 = 0.1e1 / t38 ^ 2;
t56 = t49 * t52;
t57 = t47 * t53;
t37 = t48 * t56 - t57;
t64 = t36 * t37;
t63 = t39 * t48;
t43 = t46 ^ 2;
t62 = t43 / t48 ^ 2;
t61 = t46 * t49;
t42 = 0.1e1 / (t47 ^ 2 * t62 + 0.1e1);
t60 = t47 * t42;
t54 = t37 ^ 2 * t36 + 0.1e1;
t44 = 0.1e1 / t48;
t35 = 0.1e1 / t38;
t34 = 0.1e1 / t54;
t33 = (0.1e1 + t62) * t60;
t30 = 0.1e1 / t32;
t29 = 0.1e1 / (t43 * t65 + 0.1e1);
t1 = [t44 * t42 * t61, 0, 0, t33, 0, 0; (-t30 * t59 - (-t40 * t43 * t44 * t60 + (t42 - 0.1e1) * t46 * t39) * t46 * t65) * t29, 0, 0 (t48 * t30 - (-t47 * t63 + t40 * t46 + (-t40 * t59 + t63) * t33) * t46 * t31) * t49 * t29, 0, 0; ((-t48 * t58 - t55) * t35 - (-t48 * t57 + t56) * t64) * t34, 0, 0 (-t35 * t52 + t53 * t64) * t34 * t61, t54 * t34, 0;];
Ja_rot  = t1;
