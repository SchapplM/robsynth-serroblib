% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPPRR3_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobia_rot_3_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:59
% EndTime: 2019-02-26 19:46:00
% DurationCPUTime: 0.04s
% Computational Cost: add. (69->16), mult. (180->37), div. (39->10), fcn. (276->9), ass. (0->23)
t36 = sin(pkin(6));
t40 = cos(qJ(2));
t44 = t36 * t40;
t38 = cos(pkin(6));
t39 = sin(qJ(2));
t43 = t38 * t39;
t42 = t38 * t40;
t41 = t36 ^ 2;
t37 = cos(pkin(10));
t35 = sin(pkin(10));
t34 = 0.1e1 / t40 ^ 2;
t31 = -t35 * t43 + t37 * t40;
t30 = t35 * t42 + t37 * t39;
t29 = t35 * t40 + t37 * t43;
t27 = t35 * t39 - t37 * t42;
t26 = 0.1e1 / t31 ^ 2;
t24 = atan2(-t27, -t44);
t23 = cos(t24);
t22 = sin(t24);
t21 = -t22 * t27 - t23 * t44;
t20 = 0.1e1 / t21 ^ 2;
t18 = (t29 / t40 + t39 * t27 * t34) / t36 / (0.1e1 + t27 ^ 2 / t41 * t34);
t1 = [0, t18, 0, 0, 0, 0; 0 (t31 / t21 - (t23 * t36 * t39 - t22 * t29 + (t22 * t44 - t23 * t27) * t18) * t30 * t20) / (t30 ^ 2 * t20 + 0.1e1) 0, 0, 0, 0; 0, -t30 * t35 * t36 * t26 / (t35 ^ 2 * t41 * t26 + 0.1e1) 0, 0, 0, 0;];
Ja_rot  = t1;
