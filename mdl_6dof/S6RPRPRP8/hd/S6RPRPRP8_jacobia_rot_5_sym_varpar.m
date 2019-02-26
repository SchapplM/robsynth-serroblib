% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP8_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:47:26
% EndTime: 2019-02-26 20:47:26
% DurationCPUTime: 0.07s
% Computational Cost: add. (215->20), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
t46 = qJ(3) + pkin(9);
t44 = sin(t46);
t45 = cos(t46);
t50 = cos(qJ(1));
t54 = t50 * t45;
t39 = atan2(-t54, t44);
t37 = sin(t39);
t38 = cos(t39);
t30 = -t37 * t54 + t38 * t44;
t29 = 0.1e1 / t30 ^ 2;
t48 = sin(qJ(1));
t62 = t29 * t48 ^ 2;
t47 = sin(qJ(5));
t53 = t50 * t47;
t49 = cos(qJ(5));
t56 = t48 * t49;
t36 = t44 * t56 + t53;
t34 = 0.1e1 / t36 ^ 2;
t52 = t50 * t49;
t57 = t48 * t47;
t35 = t44 * t57 - t52;
t61 = t34 * t35;
t60 = t37 * t44;
t43 = t45 ^ 2;
t59 = 0.1e1 / t44 ^ 2 * t43;
t58 = t45 * t48;
t40 = 0.1e1 / (t50 ^ 2 * t59 + 0.1e1);
t55 = t50 * t40;
t51 = t35 ^ 2 * t34 + 0.1e1;
t41 = 0.1e1 / t44;
t33 = 0.1e1 / t36;
t32 = 0.1e1 / t51;
t31 = (0.1e1 + t59) * t55;
t28 = 0.1e1 / t30;
t27 = 0.1e1 / (t43 * t62 + 0.1e1);
t1 = [t41 * t40 * t58, 0, t31, 0, 0, 0; (-t28 * t54 + (-t38 * t41 * t43 * t55 + (-t40 + 0.1e1) * t45 * t37) * t45 * t62) * t27, 0 (t44 * t28 + (t50 * t60 + t38 * t45 + (-t38 * t54 - t60) * t31) * t45 * t29) * t48 * t27, 0, 0, 0; ((t44 * t53 + t56) * t33 - (t44 * t52 - t57) * t61) * t32, 0 (t33 * t47 - t49 * t61) * t32 * t58, 0, t51 * t32, 0;];
Ja_rot  = t1;
