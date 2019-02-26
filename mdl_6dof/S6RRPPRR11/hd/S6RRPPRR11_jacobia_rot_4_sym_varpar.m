% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR11_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobia_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:19
% EndTime: 2019-02-26 21:34:19
% DurationCPUTime: 0.11s
% Computational Cost: add. (132->24), mult. (403->63), div. (67->11), fcn. (610->11), ass. (0->39)
t47 = sin(qJ(2));
t49 = cos(qJ(2));
t50 = cos(qJ(1));
t48 = sin(qJ(1));
t53 = cos(pkin(6));
t52 = t48 * t53;
t37 = t50 * t47 + t49 * t52;
t44 = sin(pkin(11));
t46 = cos(pkin(11));
t45 = sin(pkin(6));
t55 = t45 * t48;
t29 = t37 * t44 + t46 * t55;
t27 = 0.1e1 / t29 ^ 2;
t28 = -t37 * t46 + t44 * t55;
t60 = t27 * t28;
t51 = t50 * t53;
t35 = t47 * t51 + t48 * t49;
t56 = t45 * t47;
t33 = atan2(-t35, t56);
t31 = cos(t33);
t59 = t31 * t35;
t30 = sin(t33);
t24 = -t30 * t35 + t31 * t56;
t23 = 0.1e1 / t24 ^ 2;
t38 = -t47 * t52 + t50 * t49;
t58 = t38 ^ 2 * t23;
t41 = 0.1e1 / t45;
t42 = 0.1e1 / t47;
t57 = t41 * t42;
t54 = t45 * t50;
t43 = 0.1e1 / t47 ^ 2;
t34 = t48 * t47 - t49 * t51;
t32 = 0.1e1 / (0.1e1 + t35 ^ 2 / t45 ^ 2 * t43);
t26 = 0.1e1 / t29;
t25 = 0.1e1 / (t28 ^ 2 * t27 + 0.1e1);
t22 = 0.1e1 / t24;
t21 = 0.1e1 / (0.1e1 + t58);
t20 = (t35 * t43 * t49 + t34 * t42) * t41 * t32;
t1 = [-t38 * t32 * t57, t20, 0, 0, 0, 0; (-t35 * t22 - (-t30 + (t57 * t59 + t30) * t32) * t58) * t21 (-t37 * t22 - (t31 * t45 * t49 + t30 * t34 + (-t30 * t56 - t59) * t20) * t38 * t23) * t21, 0, 0, 0, 0; ((t34 * t46 + t44 * t54) * t26 - (-t34 * t44 + t46 * t54) * t60) * t25 (-t46 * t26 - t44 * t60) * t38 * t25, 0, 0, 0, 0;];
Ja_rot  = t1;
