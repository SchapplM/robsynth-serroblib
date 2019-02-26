% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR3
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
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:48
% EndTime: 2019-02-26 19:54:48
% DurationCPUTime: 0.10s
% Computational Cost: add. (206->19), mult. (316->43), div. (52->11), fcn. (458->11), ass. (0->31)
t60 = sin(pkin(11));
t61 = sin(pkin(6));
t70 = t60 * t61;
t65 = cos(qJ(2));
t69 = t61 * t65;
t63 = cos(pkin(6));
t64 = sin(qJ(2));
t68 = t63 * t64;
t67 = t63 * t65;
t62 = cos(pkin(11));
t53 = -t60 * t68 + t62 * t65;
t57 = pkin(12) + qJ(4) + qJ(5);
t55 = sin(t57);
t56 = cos(t57);
t44 = t53 * t56 + t55 * t70;
t42 = 0.1e1 / t44 ^ 2;
t43 = t53 * t55 - t56 * t70;
t66 = t43 ^ 2 * t42 + 0.1e1;
t59 = 0.1e1 / t65 ^ 2;
t52 = t60 * t67 + t62 * t64;
t51 = t60 * t65 + t62 * t68;
t49 = t60 * t64 - t62 * t67;
t47 = atan2(-t49, -t69);
t46 = cos(t47);
t45 = sin(t47);
t41 = 0.1e1 / t66;
t40 = -t45 * t49 - t46 * t69;
t39 = 0.1e1 / t40 ^ 2;
t37 = (t51 / t65 + t64 * t49 * t59) / t61 / (0.1e1 + t49 ^ 2 / t61 ^ 2 * t59);
t36 = t66 * t41;
t1 = [0, t37, 0, 0, 0, 0; 0 (t53 / t40 - (t46 * t61 * t64 - t45 * t51 + (t45 * t69 - t46 * t49) * t37) * t52 * t39) / (t52 ^ 2 * t39 + 0.1e1) 0, 0, 0, 0; 0 (-t55 / t44 + t56 * t43 * t42) * t52 * t41, 0, t36, t36, 0;];
Ja_rot  = t1;
