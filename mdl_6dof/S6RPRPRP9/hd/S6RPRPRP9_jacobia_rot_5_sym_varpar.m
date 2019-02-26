% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP9
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
% Datum: 2019-02-26 20:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP9_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:48:01
% EndTime: 2019-02-26 20:48:01
% DurationCPUTime: 0.10s
% Computational Cost: add. (134->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->37)
t44 = sin(qJ(1));
t59 = t44 ^ 2;
t43 = sin(qJ(3));
t45 = cos(qJ(3));
t46 = cos(qJ(1));
t48 = t46 * t45;
t35 = atan2(-t48, t43);
t33 = sin(t35);
t34 = cos(t35);
t27 = -t33 * t48 + t34 * t43;
t26 = 0.1e1 / t27 ^ 2;
t58 = t26 * t45;
t39 = pkin(9) + qJ(5);
t37 = sin(t39);
t50 = t46 * t37;
t38 = cos(t39);
t53 = t44 * t38;
t32 = t43 * t53 + t50;
t30 = 0.1e1 / t32 ^ 2;
t49 = t46 * t38;
t54 = t44 * t37;
t31 = t43 * t54 - t49;
t57 = t30 * t31;
t56 = t33 * t43;
t42 = t45 ^ 2;
t55 = 0.1e1 / t43 ^ 2 * t42;
t52 = t44 * t45;
t36 = 0.1e1 / (t46 ^ 2 * t55 + 0.1e1);
t51 = t46 * t36;
t47 = t31 ^ 2 * t30 + 0.1e1;
t40 = 0.1e1 / t43;
t29 = 0.1e1 / t32;
t28 = (0.1e1 + t55) * t51;
t25 = 0.1e1 / t27;
t24 = 0.1e1 / t47;
t23 = 0.1e1 / (t59 * t42 * t26 + 0.1e1);
t1 = [t40 * t36 * t52, 0, t28, 0, 0, 0; (-t25 * t48 + (-t34 * t40 * t42 * t51 + (-t36 + 0.1e1) * t45 * t33) * t59 * t58) * t23, 0 (t43 * t25 + (t46 * t56 + t34 * t45 + (-t34 * t48 - t56) * t28) * t58) * t44 * t23, 0, 0, 0; ((t43 * t50 + t53) * t29 - (t43 * t49 - t54) * t57) * t24, 0 (t29 * t37 - t38 * t57) * t24 * t52, 0, t47 * t24, 0;];
Ja_rot  = t1;
