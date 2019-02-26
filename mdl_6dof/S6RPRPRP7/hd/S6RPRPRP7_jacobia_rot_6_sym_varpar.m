% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRP7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:46:56
% EndTime: 2019-02-26 20:46:56
% DurationCPUTime: 0.10s
% Computational Cost: add. (215->20), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
t49 = qJ(3) + pkin(9);
t47 = sin(t49);
t48 = cos(t49);
t53 = cos(qJ(1));
t57 = t53 * t48;
t42 = atan2(-t57, t47);
t40 = sin(t42);
t41 = cos(t42);
t33 = -t40 * t57 + t41 * t47;
t32 = 0.1e1 / t33 ^ 2;
t51 = sin(qJ(1));
t65 = t32 * t51 ^ 2;
t50 = sin(qJ(5));
t56 = t53 * t50;
t52 = cos(qJ(5));
t59 = t51 * t52;
t39 = t47 * t59 + t56;
t37 = 0.1e1 / t39 ^ 2;
t55 = t53 * t52;
t60 = t51 * t50;
t38 = t47 * t60 - t55;
t64 = t37 * t38;
t63 = t40 * t47;
t46 = t48 ^ 2;
t62 = 0.1e1 / t47 ^ 2 * t46;
t61 = t48 * t51;
t43 = 0.1e1 / (t53 ^ 2 * t62 + 0.1e1);
t58 = t53 * t43;
t54 = t38 ^ 2 * t37 + 0.1e1;
t44 = 0.1e1 / t47;
t36 = 0.1e1 / t39;
t35 = 0.1e1 / t54;
t34 = (0.1e1 + t62) * t58;
t31 = 0.1e1 / t33;
t30 = 0.1e1 / (t46 * t65 + 0.1e1);
t1 = [t44 * t43 * t61, 0, t34, 0, 0, 0; (-t31 * t57 + (-t41 * t44 * t46 * t58 + (-t43 + 0.1e1) * t48 * t40) * t48 * t65) * t30, 0 (t47 * t31 + (t53 * t63 + t41 * t48 + (-t41 * t57 - t63) * t34) * t48 * t32) * t51 * t30, 0, 0, 0; ((t47 * t56 + t59) * t36 - (t47 * t55 - t60) * t64) * t35, 0 (t36 * t50 - t52 * t64) * t35 * t61, 0, t54 * t35, 0;];
Ja_rot  = t1;
