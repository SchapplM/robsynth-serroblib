% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP10_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_jacobia_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:48:37
% EndTime: 2019-02-26 20:48:37
% DurationCPUTime: 0.10s
% Computational Cost: add. (195->26), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->39)
t50 = sin(qJ(3));
t66 = t50 ^ 2;
t49 = sin(qJ(5));
t51 = sin(qJ(1));
t53 = cos(qJ(3));
t52 = cos(qJ(5));
t54 = cos(qJ(1));
t56 = t54 * t52;
t37 = t51 * t49 - t53 * t56;
t59 = t50 * t52;
t32 = atan2(-t37, -t59);
t30 = sin(t32);
t31 = cos(t32);
t29 = -t30 * t37 - t31 * t59;
t27 = 0.1e1 / t29 ^ 2;
t57 = t54 * t49;
t58 = t51 * t53;
t40 = t52 * t58 + t57;
t65 = t27 * t40;
t63 = t31 * t37;
t62 = t40 ^ 2 * t27;
t44 = 0.1e1 / t50;
t47 = 0.1e1 / t52;
t61 = t44 * t47;
t60 = t50 * t51;
t41 = -t49 * t58 + t56;
t36 = 0.1e1 / t41 ^ 2;
t55 = t51 ^ 2 * t66 * t36;
t48 = 0.1e1 / t52 ^ 2;
t45 = 0.1e1 / t66;
t39 = t51 * t52 + t53 * t57;
t35 = 0.1e1 / t41;
t34 = 0.1e1 / (t37 ^ 2 * t45 * t48 + 0.1e1);
t33 = 0.1e1 / (0.1e1 + t55);
t28 = (-t37 * t45 * t47 * t53 + t54) * t34;
t26 = 0.1e1 / t29;
t25 = 0.1e1 / (0.1e1 + t62);
t24 = (t37 * t48 * t49 + t39 * t47) * t44 * t34;
t1 = [t40 * t34 * t61, 0, t28, 0, t24, 0; (-t37 * t26 - (-t30 + (-t61 * t63 + t30) * t34) * t62) * t25, 0 (t28 * t63 * t65 + (-t26 * t60 - (-t31 * t53 + (t28 - t54) * t50 * t30) * t65) * t52) * t25, 0 (t41 * t26 - (t31 * t50 * t49 - t30 * t39 + (t30 * t59 - t63) * t24) * t65) * t25, 0; (-t36 * t39 * t51 - t35 * t54) * t50 * t33, 0 (-t35 * t58 + t49 * t55) * t33, 0, -t40 * t36 * t33 * t60, 0;];
Ja_rot  = t1;
