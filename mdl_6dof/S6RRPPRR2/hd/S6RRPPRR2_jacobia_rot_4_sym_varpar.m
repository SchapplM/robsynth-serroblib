% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR2_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:28:44
% EndTime: 2019-02-26 21:28:44
% DurationCPUTime: 0.07s
% Computational Cost: add. (221->21), mult. (197->53), div. (47->9), fcn. (297->9), ass. (0->35)
t38 = qJ(2) + pkin(10);
t37 = cos(t38);
t36 = sin(t38);
t41 = sin(qJ(1));
t47 = t41 * t36;
t31 = atan2(-t47, -t37);
t29 = sin(t31);
t30 = cos(t31);
t22 = -t29 * t47 - t30 * t37;
t21 = 0.1e1 / t22 ^ 2;
t42 = cos(qJ(1));
t53 = t21 * t42 ^ 2;
t40 = cos(pkin(11));
t43 = t42 * t40;
t39 = sin(pkin(11));
t46 = t41 * t39;
t28 = t37 * t43 + t46;
t26 = 0.1e1 / t28 ^ 2;
t44 = t42 * t39;
t45 = t41 * t40;
t27 = t37 * t44 - t45;
t52 = t26 * t27;
t51 = t29 * t37;
t33 = t36 ^ 2;
t50 = t33 / t37 ^ 2;
t49 = t36 * t42;
t32 = 0.1e1 / (t41 ^ 2 * t50 + 0.1e1);
t48 = t41 * t32;
t34 = 0.1e1 / t37;
t25 = 0.1e1 / t28;
t24 = 0.1e1 / (t27 ^ 2 * t26 + 0.1e1);
t23 = (0.1e1 + t50) * t48;
t20 = 0.1e1 / t22;
t19 = 0.1e1 / (t33 * t53 + 0.1e1);
t1 = [t34 * t32 * t49, t23, 0, 0, 0, 0; (-t20 * t47 - (-t30 * t33 * t34 * t48 + (t32 - 0.1e1) * t36 * t29) * t36 * t53) * t19 (t37 * t20 - (-t41 * t51 + t30 * t36 + (-t30 * t47 + t51) * t23) * t36 * t21) * t42 * t19, 0, 0, 0, 0; ((-t37 * t46 - t43) * t25 - (-t37 * t45 + t44) * t52) * t24 (-t25 * t39 + t40 * t52) * t24 * t49, 0, 0, 0, 0;];
Ja_rot  = t1;
