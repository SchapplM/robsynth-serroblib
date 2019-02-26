% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S4RPPP1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
%
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S4RPPP1_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobia_rot_3_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobia_rot_3_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:29:48
% EndTime: 2019-02-26 19:29:48
% DurationCPUTime: 0.04s
% Computational Cost: add. (58->14), mult. (143->30), div. (29->11), fcn. (226->9), ass. (0->23)
t37 = cos(pkin(4));
t36 = cos(pkin(6));
t39 = cos(qJ(1));
t40 = t39 * t36;
t34 = sin(pkin(6));
t38 = sin(qJ(1));
t41 = t38 * t34;
t25 = -t37 * t40 + t41;
t35 = sin(pkin(4));
t42 = t35 * t36;
t22 = atan2(-t25, -t42);
t20 = sin(t22);
t21 = cos(t22);
t19 = -t20 * t25 - t21 * t42;
t27 = t38 * t37 * t36 + t39 * t34;
t44 = t27 ^ 2 / t19 ^ 2;
t30 = 0.1e1 / t35;
t43 = t30 / t36;
t33 = 0.1e1 / t38 ^ 2;
t31 = 0.1e1 / t35 ^ 2;
t28 = -t37 * t41 + t40;
t23 = 0.1e1 / (0.1e1 + t25 ^ 2 * t31 / t36 ^ 2);
t1 = [t27 * t23 * t43, 0, 0, 0; (-t25 / t19 - (-t20 + (-t21 * t25 * t43 + t20) * t23) * t44) / (0.1e1 + t44) 0, 0, 0; (-t36 + (-t34 * t37 / t38 - t28 * t33) * t39) * t30 / (t28 ^ 2 * t33 * t31 + 0.1e1) 0, 0, 0;];
Ja_rot  = t1;
