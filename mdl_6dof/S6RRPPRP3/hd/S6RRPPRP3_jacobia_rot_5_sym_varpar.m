% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6RRPPRP3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_jacobia_rot_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:08
% EndTime: 2019-02-26 21:26:08
% DurationCPUTime: 0.07s
% Computational Cost: add. (87->20), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->35)
t41 = sin(qJ(2));
t42 = sin(qJ(1));
t44 = cos(qJ(2));
t50 = t42 * t44;
t35 = atan2(-t50, t41);
t33 = sin(t35);
t34 = cos(t35);
t26 = -t33 * t50 + t34 * t41;
t25 = 0.1e1 / t26 ^ 2;
t45 = cos(qJ(1));
t57 = t25 * t45 ^ 2;
t43 = cos(qJ(5));
t47 = t45 * t43;
t40 = sin(qJ(5));
t52 = t42 * t40;
t32 = t41 * t47 - t52;
t30 = 0.1e1 / t32 ^ 2;
t48 = t45 * t40;
t51 = t42 * t43;
t31 = t41 * t48 + t51;
t56 = t30 * t31;
t55 = t33 * t41;
t39 = t44 ^ 2;
t54 = 0.1e1 / t41 ^ 2 * t39;
t36 = 0.1e1 / (t42 ^ 2 * t54 + 0.1e1);
t53 = t42 * t36;
t49 = t44 * t45;
t46 = t31 ^ 2 * t30 + 0.1e1;
t37 = 0.1e1 / t41;
t29 = 0.1e1 / t32;
t28 = (0.1e1 + t54) * t53;
t27 = 0.1e1 / t46;
t24 = 0.1e1 / t26;
t23 = 0.1e1 / (t39 * t57 + 0.1e1);
t1 = [-t37 * t36 * t49, t28, 0, 0, 0, 0; (-t24 * t50 - (t34 * t37 * t39 * t53 + (t36 - 0.1e1) * t44 * t33) * t44 * t57) * t23 (-t41 * t24 - (t42 * t55 + t34 * t44 + (-t34 * t50 - t55) * t28) * t44 * t25) * t45 * t23, 0, 0, 0, 0; ((-t41 * t52 + t47) * t29 - (-t41 * t51 - t48) * t56) * t27 (t29 * t40 - t43 * t56) * t27 * t49, 0, 0, t46 * t27, 0;];
Ja_rot  = t1;
