% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPPR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:22:26
% EndTime: 2019-02-26 21:22:26
% DurationCPUTime: 0.11s
% Computational Cost: add. (264->22), mult. (224->54), div. (52->9), fcn. (332->9), ass. (0->37)
t55 = qJ(2) + pkin(9);
t51 = sin(t55);
t53 = cos(t55);
t56 = sin(qJ(1));
t61 = t56 * t53;
t45 = atan2(-t61, t51);
t43 = sin(t45);
t44 = cos(t45);
t36 = -t43 * t61 + t44 * t51;
t35 = 0.1e1 / t36 ^ 2;
t57 = cos(qJ(1));
t69 = t35 * t57 ^ 2;
t54 = pkin(10) + qJ(6);
t50 = sin(t54);
t60 = t57 * t50;
t52 = cos(t54);
t62 = t56 * t52;
t42 = t51 * t60 + t62;
t40 = 0.1e1 / t42 ^ 2;
t59 = t57 * t52;
t63 = t56 * t50;
t41 = -t51 * t59 + t63;
t68 = t40 * t41;
t67 = t43 * t51;
t49 = t53 ^ 2;
t66 = 0.1e1 / t51 ^ 2 * t49;
t65 = t53 * t57;
t46 = 0.1e1 / (t56 ^ 2 * t66 + 0.1e1);
t64 = t56 * t46;
t58 = t41 ^ 2 * t40 + 0.1e1;
t47 = 0.1e1 / t51;
t39 = 0.1e1 / t42;
t38 = (0.1e1 + t66) * t64;
t37 = 0.1e1 / t58;
t34 = 0.1e1 / t36;
t33 = 0.1e1 / (t49 * t69 + 0.1e1);
t1 = [-t47 * t46 * t65, t38, 0, 0, 0, 0; (-t34 * t61 - (t44 * t47 * t49 * t64 + (t46 - 0.1e1) * t53 * t43) * t53 * t69) * t33 (-t51 * t34 - (t56 * t67 + t44 * t53 + (-t44 * t61 - t67) * t38) * t53 * t35) * t57 * t33, 0, 0, 0, 0; ((t51 * t62 + t60) * t39 - (-t51 * t63 + t59) * t68) * t37 (-t39 * t52 - t50 * t68) * t37 * t65, 0, 0, 0, t58 * t37;];
Ja_rot  = t1;
