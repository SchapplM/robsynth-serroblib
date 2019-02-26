% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6RPRPRP3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:44:41
% EndTime: 2019-02-26 20:44:42
% DurationCPUTime: 0.10s
% Computational Cost: add. (266->22), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
t44 = qJ(1) + pkin(9);
t42 = cos(t44);
t59 = t42 ^ 2;
t49 = cos(qJ(3));
t40 = sin(t44);
t48 = sin(qJ(3));
t55 = t40 * t48;
t37 = atan2(-t55, -t49);
t35 = sin(t37);
t36 = cos(t37);
t29 = -t35 * t55 - t36 * t49;
t28 = 0.1e1 / t29 ^ 2;
t58 = t28 * t48;
t43 = pkin(10) + qJ(5);
t39 = sin(t43);
t41 = cos(t43);
t52 = t42 * t49;
t34 = t40 * t39 + t41 * t52;
t32 = 0.1e1 / t34 ^ 2;
t33 = t39 * t52 - t40 * t41;
t57 = t32 * t33;
t45 = t48 ^ 2;
t51 = t45 / t49 ^ 2;
t38 = 0.1e1 / (t40 ^ 2 * t51 + 0.1e1);
t56 = t40 * t38;
t54 = t40 * t49;
t53 = t42 * t48;
t50 = t33 ^ 2 * t32 + 0.1e1;
t46 = 0.1e1 / t49;
t31 = 0.1e1 / t34;
t30 = (0.1e1 + t51) * t56;
t27 = 0.1e1 / t29;
t26 = 0.1e1 / t50;
t25 = 0.1e1 / (t59 * t45 * t28 + 0.1e1);
t1 = [t46 * t38 * t53, 0, t30, 0, 0, 0; (-t27 * t55 - (-t36 * t45 * t46 * t56 + (t38 - 0.1e1) * t48 * t35) * t59 * t58) * t25, 0 (t49 * t27 - (-t35 * t54 + t36 * t48 + (t35 * t49 - t36 * t55) * t30) * t58) * t42 * t25, 0, 0, 0; ((-t39 * t54 - t42 * t41) * t31 - (t42 * t39 - t41 * t54) * t57) * t26, 0 (-t31 * t39 + t41 * t57) * t26 * t53, 0, t50 * t26, 0;];
Ja_rot  = t1;
