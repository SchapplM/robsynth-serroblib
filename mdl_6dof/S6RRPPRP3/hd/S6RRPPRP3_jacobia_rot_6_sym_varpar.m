% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRP3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_jacobia_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:08
% EndTime: 2019-02-26 21:26:09
% DurationCPUTime: 0.07s
% Computational Cost: add. (87->20), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->35)
t44 = sin(qJ(2));
t45 = sin(qJ(1));
t47 = cos(qJ(2));
t53 = t45 * t47;
t38 = atan2(-t53, t44);
t36 = sin(t38);
t37 = cos(t38);
t29 = -t36 * t53 + t37 * t44;
t28 = 0.1e1 / t29 ^ 2;
t48 = cos(qJ(1));
t60 = t28 * t48 ^ 2;
t46 = cos(qJ(5));
t50 = t48 * t46;
t43 = sin(qJ(5));
t55 = t45 * t43;
t35 = t44 * t50 - t55;
t33 = 0.1e1 / t35 ^ 2;
t51 = t48 * t43;
t54 = t45 * t46;
t34 = t44 * t51 + t54;
t59 = t33 * t34;
t58 = t36 * t44;
t42 = t47 ^ 2;
t57 = 0.1e1 / t44 ^ 2 * t42;
t39 = 0.1e1 / (t45 ^ 2 * t57 + 0.1e1);
t56 = t45 * t39;
t52 = t47 * t48;
t49 = t34 ^ 2 * t33 + 0.1e1;
t40 = 0.1e1 / t44;
t32 = 0.1e1 / t35;
t31 = (0.1e1 + t57) * t56;
t30 = 0.1e1 / t49;
t27 = 0.1e1 / t29;
t26 = 0.1e1 / (t42 * t60 + 0.1e1);
t1 = [-t40 * t39 * t52, t31, 0, 0, 0, 0; (-t27 * t53 - (t37 * t40 * t42 * t56 + (t39 - 0.1e1) * t47 * t36) * t47 * t60) * t26 (-t44 * t27 - (t45 * t58 + t37 * t47 + (-t37 * t53 - t58) * t31) * t47 * t28) * t48 * t26, 0, 0, 0, 0; ((-t44 * t55 + t50) * t32 - (-t44 * t54 - t51) * t59) * t30 (t32 * t43 - t46 * t59) * t30 * t52, 0, 0, t49 * t30, 0;];
Ja_rot  = t1;
