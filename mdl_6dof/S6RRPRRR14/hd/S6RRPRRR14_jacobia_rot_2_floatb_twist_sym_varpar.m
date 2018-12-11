% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S6RRPRRR14
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR14_jacobia_rot_2_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_rot_2_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_rot_2_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_2_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:20
% EndTime: 2018-12-10 18:38:20
% DurationCPUTime: 0.05s
% Computational Cost: add. (91->18), mult. (140->39), div. (26->9), fcn. (185->13), ass. (0->26)
t41 = cos(pkin(6));
t40 = sin(pkin(6));
t45 = cos(qJ(1));
t46 = t45 * t40;
t34 = atan2(t46, t41);
t29 = sin(t34);
t30 = cos(t34);
t24 = t29 * t46 + t30 * t41;
t43 = sin(qJ(1));
t49 = 0.1e1 / t24 ^ 2 * t43 ^ 2;
t38 = pkin(6) + qJ(2);
t39 = pkin(6) - qJ(2);
t31 = sin(t38) / 0.2e1 - sin(t39) / 0.2e1;
t44 = cos(qJ(2));
t28 = -t43 * t31 + t45 * t44;
t26 = 0.1e1 / t28 ^ 2;
t32 = cos(t39) / 0.2e1 + cos(t38) / 0.2e1;
t42 = sin(qJ(2));
t27 = t43 * t32 + t45 * t42;
t48 = t27 ^ 2 * t26;
t36 = t40 ^ 2;
t33 = 0.1e1 / (0.1e1 + t45 ^ 2 * t36 / t41 ^ 2);
t47 = t33 / t41;
t25 = 0.1e1 / t28;
t21 = 0.1e1 / (0.1e1 + t48);
t1 = [-t43 * t40 * t47, 0, 0, 0, 0, 0; (0.1e1 / t24 * t46 - (-t30 * t36 * t45 * t47 + (t33 - 0.1e1) * t40 * t29) * t40 * t49) / (t36 * t49 + 0.1e1) 0, 0, 0, 0, 0; ((t45 * t32 - t43 * t42) * t25 - (-t45 * t31 - t43 * t44) * t27 * t26) * t21 (t28 * t25 + t48) * t21, 0, 0, 0, 0;];
Ja_rot  = t1;
