% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6RPPRRP5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_jacobia_rot_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:32:32
% EndTime: 2019-02-26 20:32:33
% DurationCPUTime: 0.10s
% Computational Cost: add. (63->18), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->35)
t40 = sin(qJ(4));
t41 = sin(qJ(1));
t43 = cos(qJ(4));
t49 = t41 * t43;
t35 = atan2(t49, t40);
t32 = sin(t35);
t33 = cos(t35);
t25 = t32 * t49 + t33 * t40;
t24 = 0.1e1 / t25 ^ 2;
t44 = cos(qJ(1));
t56 = t24 * t44 ^ 2;
t42 = cos(qJ(5));
t46 = t44 * t42;
t39 = sin(qJ(5));
t51 = t41 * t39;
t31 = t40 * t46 - t51;
t29 = 0.1e1 / t31 ^ 2;
t47 = t44 * t39;
t50 = t41 * t42;
t30 = t40 * t47 + t50;
t55 = t29 * t30;
t54 = t32 * t40;
t38 = t43 ^ 2;
t53 = 0.1e1 / t40 ^ 2 * t38;
t34 = 0.1e1 / (t41 ^ 2 * t53 + 0.1e1);
t52 = t41 * t34;
t48 = t43 * t44;
t45 = t30 ^ 2 * t29 + 0.1e1;
t36 = 0.1e1 / t40;
t28 = 0.1e1 / t31;
t27 = (-0.1e1 - t53) * t52;
t26 = 0.1e1 / t45;
t23 = 0.1e1 / t25;
t22 = 0.1e1 / (t38 * t56 + 0.1e1);
t1 = [t36 * t34 * t48, 0, 0, t27, 0, 0; (t23 * t49 + (t33 * t36 * t38 * t52 + (-t34 + 0.1e1) * t43 * t32) * t43 * t56) * t22, 0, 0 (t40 * t23 + (-t41 * t54 + t33 * t43 + (t33 * t49 - t54) * t27) * t43 * t24) * t44 * t22, 0, 0; ((-t40 * t51 + t46) * t28 - (-t40 * t50 - t47) * t55) * t26, 0, 0 (t28 * t39 - t42 * t55) * t26 * t48, t45 * t26, 0;];
Ja_rot  = t1;
