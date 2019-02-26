% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR8_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_jacobia_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:47
% EndTime: 2019-02-26 21:41:47
% DurationCPUTime: 0.10s
% Computational Cost: add. (157->21), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->34)
t46 = cos(qJ(2));
t44 = sin(qJ(2));
t45 = sin(qJ(1));
t52 = t45 * t44;
t36 = atan2(-t52, -t46);
t34 = sin(t36);
t35 = cos(t36);
t28 = -t34 * t52 - t35 * t46;
t27 = 0.1e1 / t28 ^ 2;
t47 = cos(qJ(1));
t57 = t27 * t47 ^ 2;
t40 = pkin(10) + qJ(4);
t38 = sin(t40);
t39 = cos(t40);
t49 = t47 * t39;
t33 = t45 * t38 + t46 * t49;
t31 = 0.1e1 / t33 ^ 2;
t50 = t47 * t38;
t32 = -t45 * t39 + t46 * t50;
t56 = t31 * t32;
t41 = t44 ^ 2;
t55 = t41 / t46 ^ 2;
t54 = t44 * t47;
t37 = 0.1e1 / (t45 ^ 2 * t55 + 0.1e1);
t53 = t45 * t37;
t51 = t45 * t46;
t48 = t32 ^ 2 * t31 + 0.1e1;
t42 = 0.1e1 / t46;
t30 = 0.1e1 / t33;
t29 = (0.1e1 + t55) * t53;
t26 = 0.1e1 / t28;
t25 = 0.1e1 / t48;
t24 = 0.1e1 / (t41 * t57 + 0.1e1);
t1 = [t42 * t37 * t54, t29, 0, 0, 0, 0; (-t26 * t52 - (-t35 * t41 * t42 * t53 + (t37 - 0.1e1) * t44 * t34) * t44 * t57) * t24 (t46 * t26 - (-t34 * t51 + t35 * t44 + (t34 * t46 - t35 * t52) * t29) * t44 * t27) * t47 * t24, 0, 0, 0, 0; ((-t38 * t51 - t49) * t30 - (-t39 * t51 + t50) * t56) * t25 (-t30 * t38 + t39 * t56) * t25 * t54, 0, t48 * t25, 0, 0;];
Ja_rot  = t1;
