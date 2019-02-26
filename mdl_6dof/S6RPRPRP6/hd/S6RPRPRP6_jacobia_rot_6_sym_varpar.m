% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:46:27
% EndTime: 2019-02-26 20:46:27
% DurationCPUTime: 0.07s
% Computational Cost: add. (216->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
t48 = pkin(9) + qJ(3);
t46 = sin(t48);
t47 = cos(t48);
t50 = sin(qJ(1));
t58 = t50 * t47;
t41 = atan2(-t58, t46);
t39 = sin(t41);
t40 = cos(t41);
t32 = -t39 * t58 + t40 * t46;
t31 = 0.1e1 / t32 ^ 2;
t52 = cos(qJ(1));
t64 = t31 * t52 ^ 2;
t49 = sin(qJ(5));
t55 = t52 * t49;
t51 = cos(qJ(5));
t56 = t50 * t51;
t38 = t46 * t55 + t56;
t36 = 0.1e1 / t38 ^ 2;
t54 = t52 * t51;
t57 = t50 * t49;
t37 = -t46 * t54 + t57;
t63 = t36 * t37;
t62 = t39 * t46;
t45 = t47 ^ 2;
t61 = 0.1e1 / t46 ^ 2 * t45;
t60 = t47 * t52;
t42 = 0.1e1 / (t50 ^ 2 * t61 + 0.1e1);
t59 = t50 * t42;
t53 = t37 ^ 2 * t36 + 0.1e1;
t43 = 0.1e1 / t46;
t35 = 0.1e1 / t38;
t34 = 0.1e1 / t53;
t33 = (0.1e1 + t61) * t59;
t30 = 0.1e1 / t32;
t29 = 0.1e1 / (t45 * t64 + 0.1e1);
t1 = [-t43 * t42 * t60, 0, t33, 0, 0, 0; (-t30 * t58 - (t40 * t43 * t45 * t59 + (t42 - 0.1e1) * t47 * t39) * t47 * t64) * t29, 0 (-t46 * t30 - (t50 * t62 + t40 * t47 + (-t40 * t58 - t62) * t33) * t47 * t31) * t52 * t29, 0, 0, 0; ((t46 * t56 + t55) * t35 - (-t46 * t57 + t54) * t63) * t34, 0 (-t35 * t51 - t49 * t63) * t34 * t60, 0, t53 * t34, 0;];
Ja_rot  = t1;
