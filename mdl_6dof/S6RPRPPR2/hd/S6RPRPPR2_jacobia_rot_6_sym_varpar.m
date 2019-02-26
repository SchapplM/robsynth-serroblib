% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPPR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:39:46
% EndTime: 2019-02-26 20:39:46
% DurationCPUTime: 0.10s
% Computational Cost: add. (325->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
t50 = qJ(3) + pkin(10);
t46 = sin(t50);
t51 = qJ(1) + pkin(9);
t47 = sin(t51);
t48 = cos(t50);
t60 = t47 * t48;
t41 = atan2(-t60, t46);
t39 = sin(t41);
t40 = cos(t41);
t33 = -t39 * t60 + t40 * t46;
t31 = 0.1e1 / t33 ^ 2;
t49 = cos(t51);
t65 = t31 * t49 ^ 2;
t52 = sin(qJ(6));
t56 = t49 * t52;
t53 = cos(qJ(6));
t58 = t47 * t53;
t38 = t46 * t56 + t58;
t36 = 0.1e1 / t38 ^ 2;
t55 = t49 * t53;
t59 = t47 * t52;
t37 = -t46 * t55 + t59;
t64 = t36 * t37;
t63 = t39 * t46;
t45 = t48 ^ 2;
t62 = 0.1e1 / t46 ^ 2 * t45;
t42 = 0.1e1 / (t47 ^ 2 * t62 + 0.1e1);
t61 = t47 * t42;
t57 = t48 * t49;
t54 = t37 ^ 2 * t36 + 0.1e1;
t43 = 0.1e1 / t46;
t35 = 0.1e1 / t38;
t34 = 0.1e1 / t54;
t32 = (0.1e1 + t62) * t61;
t30 = 0.1e1 / t33;
t29 = 0.1e1 / (t45 * t65 + 0.1e1);
t1 = [-t43 * t42 * t57, 0, t32, 0, 0, 0; (-t30 * t60 - (t40 * t43 * t45 * t61 + (t42 - 0.1e1) * t48 * t39) * t48 * t65) * t29, 0 (-t46 * t30 - (t47 * t63 + t40 * t48 + (-t40 * t60 - t63) * t32) * t48 * t31) * t49 * t29, 0, 0, 0; ((t46 * t58 + t56) * t35 - (-t46 * t59 + t55) * t64) * t34, 0 (-t35 * t53 - t52 * t64) * t34 * t57, 0, 0, t54 * t34;];
Ja_rot  = t1;
