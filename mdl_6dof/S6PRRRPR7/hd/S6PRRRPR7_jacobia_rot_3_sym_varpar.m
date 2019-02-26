% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRRRPR7
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR7_jacobia_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobia_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobia_rot_3_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:04
% EndTime: 2019-02-26 20:14:04
% DurationCPUTime: 0.11s
% Computational Cost: add. (179->22), mult. (525->56), div. (35->9), fcn. (736->13), ass. (0->37)
t45 = sin(pkin(12));
t48 = cos(pkin(12));
t54 = cos(qJ(2));
t50 = cos(pkin(6));
t52 = sin(qJ(2));
t58 = t50 * t52;
t43 = -t45 * t58 + t48 * t54;
t51 = sin(qJ(3));
t63 = t43 * t51;
t53 = cos(qJ(3));
t62 = t43 * t53;
t46 = sin(pkin(7));
t47 = sin(pkin(6));
t61 = t46 * t47;
t49 = cos(pkin(7));
t60 = t47 * t49;
t59 = t47 * t52;
t57 = t50 * t54;
t42 = -t45 * t57 - t48 * t52;
t55 = t42 * t49 + t45 * t61;
t32 = t51 * t55 + t62;
t30 = 0.1e1 / t32 ^ 2;
t31 = -t53 * t55 + t63;
t56 = t31 ^ 2 * t30 + 0.1e1;
t41 = -t45 * t54 - t48 * t58;
t40 = t50 * t49 - t54 * t61;
t39 = 0.1e1 / t40 ^ 2;
t38 = -t42 * t46 + t45 * t60;
t37 = (-t45 * t52 + t48 * t57) * t46 + t48 * t60;
t36 = atan2(t37, t40);
t34 = cos(t36);
t33 = sin(t36);
t29 = 0.1e1 / t56;
t28 = t33 * t37 + t34 * t40;
t27 = 0.1e1 / t28 ^ 2;
t25 = (t41 / t40 - t37 * t39 * t59) * t46 / (t37 ^ 2 * t39 + 0.1e1);
t1 = [0, t25, 0, 0, 0, 0; 0 (t43 * t46 / t28 - ((t33 * t41 + t34 * t59) * t46 + (-t33 * t40 + t34 * t37) * t25) * t38 * t27) / (t38 ^ 2 * t27 + 0.1e1) 0, 0, 0, 0; 0 ((t42 * t51 + t49 * t62) / t32 - (t42 * t53 - t49 * t63) * t31 * t30) * t29, t56 * t29, 0, 0, 0;];
Ja_rot  = t1;
