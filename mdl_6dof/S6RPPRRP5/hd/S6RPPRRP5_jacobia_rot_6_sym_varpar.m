% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRP5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_jacobia_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:32:32
% EndTime: 2019-02-26 20:32:33
% DurationCPUTime: 0.10s
% Computational Cost: add. (63->18), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->35)
t43 = sin(qJ(4));
t44 = sin(qJ(1));
t46 = cos(qJ(4));
t52 = t44 * t46;
t38 = atan2(t52, t43);
t35 = sin(t38);
t36 = cos(t38);
t28 = t35 * t52 + t36 * t43;
t27 = 0.1e1 / t28 ^ 2;
t47 = cos(qJ(1));
t59 = t27 * t47 ^ 2;
t45 = cos(qJ(5));
t49 = t47 * t45;
t42 = sin(qJ(5));
t54 = t44 * t42;
t34 = t43 * t49 - t54;
t32 = 0.1e1 / t34 ^ 2;
t50 = t47 * t42;
t53 = t44 * t45;
t33 = t43 * t50 + t53;
t58 = t32 * t33;
t57 = t35 * t43;
t41 = t46 ^ 2;
t56 = 0.1e1 / t43 ^ 2 * t41;
t37 = 0.1e1 / (t44 ^ 2 * t56 + 0.1e1);
t55 = t44 * t37;
t51 = t46 * t47;
t48 = t33 ^ 2 * t32 + 0.1e1;
t39 = 0.1e1 / t43;
t31 = 0.1e1 / t34;
t30 = (-0.1e1 - t56) * t55;
t29 = 0.1e1 / t48;
t26 = 0.1e1 / t28;
t25 = 0.1e1 / (t41 * t59 + 0.1e1);
t1 = [t39 * t37 * t51, 0, 0, t30, 0, 0; (t26 * t52 + (t36 * t39 * t41 * t55 + (-t37 + 0.1e1) * t46 * t35) * t46 * t59) * t25, 0, 0 (t43 * t26 + (-t44 * t57 + t36 * t46 + (t36 * t52 - t57) * t30) * t46 * t27) * t47 * t25, 0, 0; ((-t43 * t54 + t49) * t31 - (-t43 * t53 - t50) * t58) * t29, 0, 0 (t31 * t42 - t45 * t58) * t29 * t51, t48 * t29, 0;];
Ja_rot  = t1;
