% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRP2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:25:21
% EndTime: 2019-02-26 21:25:21
% DurationCPUTime: 0.10s
% Computational Cost: add. (216->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
t47 = qJ(2) + pkin(9);
t45 = sin(t47);
t46 = cos(t47);
t49 = sin(qJ(1));
t57 = t49 * t46;
t40 = atan2(-t57, t45);
t38 = sin(t40);
t39 = cos(t40);
t31 = -t38 * t57 + t39 * t45;
t30 = 0.1e1 / t31 ^ 2;
t51 = cos(qJ(1));
t63 = t30 * t51 ^ 2;
t48 = sin(qJ(5));
t54 = t51 * t48;
t50 = cos(qJ(5));
t55 = t49 * t50;
t37 = t45 * t54 + t55;
t35 = 0.1e1 / t37 ^ 2;
t53 = t51 * t50;
t56 = t49 * t48;
t36 = -t45 * t53 + t56;
t62 = t35 * t36;
t61 = t38 * t45;
t44 = t46 ^ 2;
t60 = 0.1e1 / t45 ^ 2 * t44;
t59 = t46 * t51;
t41 = 0.1e1 / (t49 ^ 2 * t60 + 0.1e1);
t58 = t49 * t41;
t52 = t36 ^ 2 * t35 + 0.1e1;
t42 = 0.1e1 / t45;
t34 = 0.1e1 / t37;
t33 = 0.1e1 / t52;
t32 = (0.1e1 + t60) * t58;
t29 = 0.1e1 / t31;
t28 = 0.1e1 / (t44 * t63 + 0.1e1);
t1 = [-t42 * t41 * t59, t32, 0, 0, 0, 0; (-t29 * t57 - (t39 * t42 * t44 * t58 + (t41 - 0.1e1) * t46 * t38) * t46 * t63) * t28 (-t45 * t29 - (t49 * t61 + t39 * t46 + (-t39 * t57 - t61) * t32) * t46 * t30) * t51 * t28, 0, 0, 0, 0; ((t45 * t55 + t54) * t34 - (-t45 * t56 + t53) * t62) * t33 (-t34 * t50 - t48 * t62) * t33 * t59, 0, 0, t52 * t33, 0;];
Ja_rot  = t1;
