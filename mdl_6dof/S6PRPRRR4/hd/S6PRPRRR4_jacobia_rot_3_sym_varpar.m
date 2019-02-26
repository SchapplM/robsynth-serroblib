% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRPRRR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR4_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_jacobia_rot_3_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:55:23
% EndTime: 2019-02-26 19:55:23
% DurationCPUTime: 0.09s
% Computational Cost: add. (84->18), mult. (220->41), div. (42->11), fcn. (332->11), ass. (0->27)
t39 = sin(pkin(11));
t40 = sin(pkin(6));
t49 = t39 * t40;
t45 = cos(qJ(2));
t48 = t40 * t45;
t43 = cos(pkin(6));
t44 = sin(qJ(2));
t47 = t43 * t44;
t46 = t43 * t45;
t42 = cos(pkin(11));
t41 = cos(pkin(12));
t38 = sin(pkin(12));
t37 = 0.1e1 / t45 ^ 2;
t34 = -t39 * t47 + t42 * t45;
t33 = t39 * t46 + t42 * t44;
t32 = t39 * t45 + t42 * t47;
t30 = t39 * t44 - t42 * t46;
t28 = atan2(-t30, -t48);
t27 = cos(t28);
t26 = sin(t28);
t25 = t34 * t41 + t38 * t49;
t24 = t34 * t38 - t41 * t49;
t23 = 0.1e1 / t25 ^ 2;
t21 = -t26 * t30 - t27 * t48;
t20 = 0.1e1 / t21 ^ 2;
t18 = (t32 / t45 + t44 * t30 * t37) / t40 / (0.1e1 + t30 ^ 2 / t40 ^ 2 * t37);
t1 = [0, t18, 0, 0, 0, 0; 0 (t34 / t21 - (t27 * t40 * t44 - t26 * t32 + (t26 * t48 - t27 * t30) * t18) * t33 * t20) / (t33 ^ 2 * t20 + 0.1e1) 0, 0, 0, 0; 0 (-t38 / t25 + t41 * t24 * t23) * t33 / (t24 ^ 2 * t23 + 0.1e1) 0, 0, 0, 0;];
Ja_rot  = t1;
