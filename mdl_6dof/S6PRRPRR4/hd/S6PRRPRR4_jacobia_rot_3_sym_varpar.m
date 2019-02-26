% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR4_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobia_rot_3_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:52
% EndTime: 2019-02-26 20:05:52
% DurationCPUTime: 0.06s
% Computational Cost: add. (101->18), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->29)
t40 = sin(pkin(11));
t41 = sin(pkin(6));
t52 = t40 * t41;
t47 = cos(qJ(2));
t51 = t41 * t47;
t43 = cos(pkin(6));
t45 = sin(qJ(2));
t50 = t43 * t45;
t49 = t43 * t47;
t42 = cos(pkin(11));
t36 = -t40 * t50 + t42 * t47;
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t27 = t36 * t46 + t44 * t52;
t25 = 0.1e1 / t27 ^ 2;
t26 = t36 * t44 - t46 * t52;
t48 = t26 ^ 2 * t25 + 0.1e1;
t39 = 0.1e1 / t47 ^ 2;
t35 = t40 * t49 + t42 * t45;
t34 = t40 * t47 + t42 * t50;
t32 = t40 * t45 - t42 * t49;
t30 = atan2(-t32, -t51);
t29 = cos(t30);
t28 = sin(t30);
t24 = 0.1e1 / t48;
t23 = -t28 * t32 - t29 * t51;
t22 = 0.1e1 / t23 ^ 2;
t20 = (t34 / t47 + t45 * t32 * t39) / t41 / (0.1e1 + t32 ^ 2 / t41 ^ 2 * t39);
t1 = [0, t20, 0, 0, 0, 0; 0 (t36 / t23 - (t29 * t41 * t45 - t28 * t34 + (t28 * t51 - t29 * t32) * t20) * t35 * t22) / (t35 ^ 2 * t22 + 0.1e1) 0, 0, 0, 0; 0 (-t44 / t27 + t46 * t26 * t25) * t35 * t24, t48 * t24, 0, 0, 0;];
Ja_rot  = t1;
