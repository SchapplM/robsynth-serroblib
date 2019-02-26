% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP5_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_jacobia_rot_3_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:20
% EndTime: 2019-02-26 21:48:20
% DurationCPUTime: 0.06s
% Computational Cost: add. (73->16), mult. (212->41), div. (26->9), fcn. (311->11), ass. (0->28)
t48 = cos(pkin(6));
t45 = sin(pkin(11));
t47 = cos(pkin(11));
t49 = sin(qJ(2));
t51 = cos(qJ(2));
t54 = t51 * t45 + t49 * t47;
t35 = t54 * t48;
t36 = t49 * t45 - t51 * t47;
t50 = sin(qJ(1));
t52 = cos(qJ(1));
t31 = -t50 * t35 - t52 * t36;
t28 = 0.1e1 / t31 ^ 2;
t53 = t36 * t48;
t29 = t50 * t53 - t52 * t54;
t58 = t29 ^ 2 * t28;
t46 = sin(pkin(6));
t55 = t52 * t46;
t41 = atan2(t55, t48);
t38 = sin(t41);
t39 = cos(t41);
t33 = t38 * t55 + t39 * t48;
t57 = 0.1e1 / t33 ^ 2 * t50 ^ 2;
t43 = t46 ^ 2;
t40 = 0.1e1 / (0.1e1 + t52 ^ 2 * t43 / t48 ^ 2);
t56 = t40 / t48;
t27 = 0.1e1 / t31;
t25 = 0.1e1 / (0.1e1 + t58);
t1 = [-t50 * t46 * t56, 0, 0, 0, 0, 0; (0.1e1 / t33 * t55 - (-t39 * t43 * t52 * t56 + (t40 - 0.1e1) * t46 * t38) * t46 * t57) / (t43 * t57 + 0.1e1) 0, 0, 0, 0, 0; ((-t50 * t54 - t52 * t53) * t27 + (-t52 * t35 + t50 * t36) * t29 * t28) * t25 (t31 * t27 + t58) * t25, 0, 0, 0, 0;];
Ja_rot  = t1;
