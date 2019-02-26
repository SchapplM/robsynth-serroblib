% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPPRR1_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobia_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:35
% EndTime: 2019-02-26 19:44:35
% DurationCPUTime: 0.08s
% Computational Cost: add. (191->18), mult. (543->44), div. (30->9), fcn. (768->13), ass. (0->32)
t54 = sin(pkin(10));
t55 = sin(pkin(6));
t64 = t54 * t55;
t59 = cos(pkin(6));
t53 = sin(pkin(11));
t57 = cos(pkin(11));
t60 = sin(qJ(2));
t61 = cos(qJ(2));
t63 = t61 * t53 + t60 * t57;
t49 = t63 * t59;
t50 = t60 * t53 - t61 * t57;
t58 = cos(pkin(10));
t44 = -t54 * t49 - t58 * t50;
t62 = t50 * t59;
t56 = cos(pkin(12));
t52 = sin(pkin(12));
t48 = t63 * t55;
t47 = t50 * t55;
t46 = 0.1e1 / t47 ^ 2;
t42 = t54 * t62 - t58 * t63;
t41 = -t54 * t63 - t58 * t62;
t40 = -t58 * t49 + t54 * t50;
t39 = t44 * t56 + t52 * t64;
t38 = t44 * t52 - t56 * t64;
t37 = 0.1e1 / t39 ^ 2;
t36 = atan2(t41, t47);
t34 = cos(t36);
t33 = sin(t36);
t31 = t33 * t41 + t34 * t47;
t30 = 0.1e1 / t31 ^ 2;
t28 = (t40 / t47 - t48 * t41 * t46) / (t41 ^ 2 * t46 + 0.1e1);
t1 = [0, t28, 0, 0, 0, 0; 0 (t44 / t31 + (t33 * t40 + t34 * t48 + (-t33 * t47 + t34 * t41) * t28) * t42 * t30) / (t42 ^ 2 * t30 + 0.1e1) 0, 0, 0, 0; 0 (t52 / t39 - t56 * t38 * t37) * t42 / (t38 ^ 2 * t37 + 0.1e1) 0, 0, 0, 0;];
Ja_rot  = t1;
