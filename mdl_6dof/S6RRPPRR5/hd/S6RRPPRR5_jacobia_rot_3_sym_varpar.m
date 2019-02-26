% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR5_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobia_rot_3_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:57
% EndTime: 2019-02-26 21:30:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (128->20), mult. (330->51), div. (64->11), fcn. (506->9), ass. (0->33)
t40 = sin(qJ(2));
t41 = sin(qJ(1));
t42 = cos(qJ(2));
t43 = cos(qJ(1));
t47 = cos(pkin(6));
t45 = t43 * t47;
t30 = t41 * t40 - t42 * t45;
t39 = sin(pkin(6));
t48 = t39 * t42;
t26 = atan2(-t30, -t48);
t25 = cos(t26);
t53 = t25 * t30;
t46 = t41 * t47;
t34 = -t40 * t46 + t43 * t42;
t29 = 0.1e1 / t34 ^ 2;
t44 = t39 ^ 2;
t52 = 0.1e1 / (t41 ^ 2 * t44 * t29 + 0.1e1) * t39;
t51 = t29 * t41;
t24 = sin(t26);
t23 = -t24 * t30 - t25 * t48;
t22 = 0.1e1 / t23 ^ 2;
t33 = t43 * t40 + t42 * t46;
t50 = t33 ^ 2 * t22;
t36 = 0.1e1 / t39;
t37 = 0.1e1 / t42;
t49 = t36 * t37;
t38 = 0.1e1 / t42 ^ 2;
t32 = t40 * t45 + t41 * t42;
t28 = 0.1e1 / (0.1e1 + t30 ^ 2 / t44 * t38);
t21 = 0.1e1 / t23;
t20 = 0.1e1 / (0.1e1 + t50);
t19 = (t30 * t38 * t40 + t32 * t37) * t36 * t28;
t1 = [t33 * t28 * t49, t19, 0, 0, 0, 0; (-t30 * t21 - (-t24 + (-t49 * t53 + t24) * t28) * t50) * t20 (t34 * t21 - (t25 * t39 * t40 - t24 * t32 + (t24 * t48 - t53) * t19) * t33 * t22) * t20, 0, 0, 0, 0; (-t43 / t34 - t32 * t51) * t52, -t33 * t51 * t52, 0, 0, 0, 0;];
Ja_rot  = t1;
