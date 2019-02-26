% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP5
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
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:45:52
% EndTime: 2019-02-26 20:45:53
% DurationCPUTime: 0.10s
% Computational Cost: add. (286->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
t49 = pkin(9) + qJ(3);
t47 = cos(t49);
t45 = sin(t49);
t50 = sin(qJ(1));
t56 = t50 * t45;
t39 = atan2(-t56, -t47);
t37 = sin(t39);
t38 = cos(t39);
t30 = -t37 * t56 - t38 * t47;
t29 = 0.1e1 / t30 ^ 2;
t51 = cos(qJ(1));
t63 = t29 * t51 ^ 2;
t48 = pkin(10) + qJ(5);
t46 = cos(t48);
t53 = t51 * t46;
t44 = sin(t48);
t57 = t50 * t44;
t36 = t47 * t53 + t57;
t34 = 0.1e1 / t36 ^ 2;
t54 = t51 * t44;
t55 = t50 * t46;
t35 = t47 * t54 - t55;
t62 = t34 * t35;
t61 = t37 * t47;
t41 = t45 ^ 2;
t60 = t41 / t47 ^ 2;
t59 = t45 * t51;
t40 = 0.1e1 / (t50 ^ 2 * t60 + 0.1e1);
t58 = t50 * t40;
t52 = t35 ^ 2 * t34 + 0.1e1;
t42 = 0.1e1 / t47;
t33 = 0.1e1 / t36;
t32 = (0.1e1 + t60) * t58;
t31 = 0.1e1 / t52;
t28 = 0.1e1 / t30;
t27 = 0.1e1 / (t41 * t63 + 0.1e1);
t1 = [t42 * t40 * t59, 0, t32, 0, 0, 0; (-t28 * t56 - (-t38 * t41 * t42 * t58 + (t40 - 0.1e1) * t45 * t37) * t45 * t63) * t27, 0 (t47 * t28 - (-t50 * t61 + t38 * t45 + (-t38 * t56 + t61) * t32) * t45 * t29) * t51 * t27, 0, 0, 0; ((-t47 * t57 - t53) * t33 - (-t47 * t55 + t54) * t62) * t31, 0 (-t33 * t44 + t46 * t62) * t31 * t59, 0, t52 * t31, 0;];
Ja_rot  = t1;
