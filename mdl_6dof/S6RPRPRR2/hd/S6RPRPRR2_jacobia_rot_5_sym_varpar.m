% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:49:48
% EndTime: 2019-02-26 20:49:48
% DurationCPUTime: 0.10s
% Computational Cost: add. (347->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
t49 = qJ(3) + pkin(11);
t47 = cos(t49);
t45 = sin(t49);
t50 = qJ(1) + pkin(10);
t46 = sin(t50);
t58 = t46 * t45;
t40 = atan2(-t58, -t47);
t38 = sin(t40);
t39 = cos(t40);
t31 = -t38 * t58 - t39 * t47;
t30 = 0.1e1 / t31 ^ 2;
t48 = cos(t50);
t64 = t30 * t48 ^ 2;
t52 = cos(qJ(5));
t54 = t48 * t52;
t51 = sin(qJ(5));
t57 = t46 * t51;
t37 = t47 * t54 + t57;
t35 = 0.1e1 / t37 ^ 2;
t55 = t48 * t51;
t56 = t46 * t52;
t36 = t47 * t55 - t56;
t63 = t35 * t36;
t62 = t38 * t47;
t42 = t45 ^ 2;
t61 = t42 / t47 ^ 2;
t60 = t45 * t48;
t41 = 0.1e1 / (t46 ^ 2 * t61 + 0.1e1);
t59 = t46 * t41;
t53 = t36 ^ 2 * t35 + 0.1e1;
t43 = 0.1e1 / t47;
t34 = 0.1e1 / t37;
t33 = 0.1e1 / t53;
t32 = (0.1e1 + t61) * t59;
t29 = 0.1e1 / t31;
t28 = 0.1e1 / (t42 * t64 + 0.1e1);
t1 = [t43 * t41 * t60, 0, t32, 0, 0, 0; (-t29 * t58 - (-t39 * t42 * t43 * t59 + (t41 - 0.1e1) * t45 * t38) * t45 * t64) * t28, 0 (t47 * t29 - (-t46 * t62 + t39 * t45 + (-t39 * t58 + t62) * t32) * t45 * t30) * t48 * t28, 0, 0, 0; ((-t47 * t57 - t54) * t34 - (-t47 * t56 + t55) * t63) * t33, 0 (-t34 * t51 + t52 * t63) * t33 * t60, 0, t53 * t33, 0;];
Ja_rot  = t1;
