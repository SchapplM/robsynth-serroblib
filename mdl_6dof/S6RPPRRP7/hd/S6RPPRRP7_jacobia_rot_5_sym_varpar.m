% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRP7_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:33:36
% EndTime: 2019-02-26 20:33:36
% DurationCPUTime: 0.10s
% Computational Cost: add. (215->20), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
t43 = pkin(9) + qJ(4);
t41 = sin(t43);
t42 = cos(t43);
t47 = cos(qJ(1));
t51 = t47 * t42;
t36 = atan2(-t51, t41);
t34 = sin(t36);
t35 = cos(t36);
t27 = -t34 * t51 + t35 * t41;
t26 = 0.1e1 / t27 ^ 2;
t45 = sin(qJ(1));
t59 = t26 * t45 ^ 2;
t44 = sin(qJ(5));
t50 = t47 * t44;
t46 = cos(qJ(5));
t53 = t45 * t46;
t33 = t41 * t53 + t50;
t31 = 0.1e1 / t33 ^ 2;
t49 = t47 * t46;
t54 = t45 * t44;
t32 = t41 * t54 - t49;
t58 = t31 * t32;
t57 = t34 * t41;
t40 = t42 ^ 2;
t56 = 0.1e1 / t41 ^ 2 * t40;
t55 = t42 * t45;
t37 = 0.1e1 / (t47 ^ 2 * t56 + 0.1e1);
t52 = t47 * t37;
t48 = t32 ^ 2 * t31 + 0.1e1;
t38 = 0.1e1 / t41;
t30 = 0.1e1 / t33;
t29 = 0.1e1 / t48;
t28 = (0.1e1 + t56) * t52;
t25 = 0.1e1 / t27;
t24 = 0.1e1 / (t40 * t59 + 0.1e1);
t1 = [t38 * t37 * t55, 0, 0, t28, 0, 0; (-t25 * t51 + (-t35 * t38 * t40 * t52 + (-t37 + 0.1e1) * t42 * t34) * t42 * t59) * t24, 0, 0 (t41 * t25 + (t47 * t57 + t35 * t42 + (-t35 * t51 - t57) * t28) * t42 * t26) * t45 * t24, 0, 0; ((t41 * t50 + t53) * t30 - (t41 * t49 - t54) * t58) * t29, 0, 0 (t30 * t44 - t46 * t58) * t29 * t55, t48 * t29, 0;];
Ja_rot  = t1;
