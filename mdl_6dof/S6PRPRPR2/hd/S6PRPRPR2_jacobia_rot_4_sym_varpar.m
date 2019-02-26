% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR2_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobia_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:03
% EndTime: 2019-02-26 19:47:03
% DurationCPUTime: 0.08s
% Computational Cost: add. (222->18), mult. (626->46), div. (35->9), fcn. (880->13), ass. (0->34)
t56 = sin(pkin(10));
t57 = sin(pkin(6));
t68 = t56 * t57;
t60 = cos(pkin(6));
t55 = sin(pkin(11));
t58 = cos(pkin(11));
t62 = sin(qJ(2));
t64 = cos(qJ(2));
t66 = t64 * t55 + t62 * t58;
t52 = t66 * t60;
t53 = t62 * t55 - t64 * t58;
t59 = cos(pkin(10));
t47 = -t56 * t52 - t59 * t53;
t61 = sin(qJ(4));
t63 = cos(qJ(4));
t42 = t47 * t63 + t61 * t68;
t40 = 0.1e1 / t42 ^ 2;
t41 = t47 * t61 - t63 * t68;
t67 = t41 ^ 2 * t40 + 0.1e1;
t65 = t53 * t60;
t51 = t66 * t57;
t50 = t53 * t57;
t49 = 0.1e1 / t50 ^ 2;
t45 = t56 * t65 - t59 * t66;
t44 = -t56 * t66 - t59 * t65;
t43 = -t59 * t52 + t56 * t53;
t39 = atan2(t44, t50);
t37 = cos(t39);
t36 = sin(t39);
t35 = 0.1e1 / t67;
t34 = t36 * t44 + t37 * t50;
t33 = 0.1e1 / t34 ^ 2;
t31 = (t43 / t50 - t51 * t44 * t49) / (t44 ^ 2 * t49 + 0.1e1);
t1 = [0, t31, 0, 0, 0, 0; 0 (t47 / t34 + (t36 * t43 + t37 * t51 + (-t36 * t50 + t37 * t44) * t31) * t45 * t33) / (t45 ^ 2 * t33 + 0.1e1) 0, 0, 0, 0; 0 (t61 / t42 - t63 * t41 * t40) * t45 * t35, 0, t67 * t35, 0, 0;];
Ja_rot  = t1;
