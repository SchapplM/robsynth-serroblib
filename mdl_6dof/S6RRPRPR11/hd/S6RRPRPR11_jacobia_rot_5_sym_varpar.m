% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR11
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR11_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:36
% EndTime: 2019-02-26 21:43:36
% DurationCPUTime: 0.07s
% Computational Cost: add. (135->21), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->36)
t47 = sin(qJ(2));
t48 = sin(qJ(1));
t49 = cos(qJ(2));
t55 = t48 * t49;
t39 = atan2(-t55, t47);
t37 = sin(t39);
t38 = cos(t39);
t31 = -t37 * t55 + t38 * t47;
t30 = 0.1e1 / t31 ^ 2;
t50 = cos(qJ(1));
t62 = t30 * t50 ^ 2;
t43 = qJ(4) + pkin(10);
t41 = sin(t43);
t53 = t50 * t41;
t42 = cos(t43);
t56 = t48 * t42;
t36 = t47 * t53 + t56;
t34 = 0.1e1 / t36 ^ 2;
t52 = t50 * t42;
t57 = t48 * t41;
t35 = -t47 * t52 + t57;
t61 = t34 * t35;
t60 = t37 * t47;
t46 = t49 ^ 2;
t59 = 0.1e1 / t47 ^ 2 * t46;
t40 = 0.1e1 / (t48 ^ 2 * t59 + 0.1e1);
t58 = t48 * t40;
t54 = t49 * t50;
t51 = t35 ^ 2 * t34 + 0.1e1;
t44 = 0.1e1 / t47;
t33 = 0.1e1 / t36;
t32 = (0.1e1 + t59) * t58;
t29 = 0.1e1 / t31;
t28 = 0.1e1 / t51;
t27 = 0.1e1 / (t46 * t62 + 0.1e1);
t1 = [-t44 * t40 * t54, t32, 0, 0, 0, 0; (-t29 * t55 - (t38 * t44 * t46 * t58 + (t40 - 0.1e1) * t49 * t37) * t49 * t62) * t27 (-t47 * t29 - (t48 * t60 + t38 * t49 + (-t38 * t55 - t60) * t32) * t49 * t30) * t50 * t27, 0, 0, 0, 0; ((t47 * t56 + t53) * t33 - (-t47 * t57 + t52) * t61) * t28 (-t33 * t42 - t41 * t61) * t28 * t54, 0, t51 * t28, 0, 0;];
Ja_rot  = t1;
