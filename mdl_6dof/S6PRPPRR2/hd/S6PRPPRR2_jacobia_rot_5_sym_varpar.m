% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPPRR2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:12
% EndTime: 2019-02-26 19:45:13
% DurationCPUTime: 0.11s
% Computational Cost: add. (222->19), mult. (626->46), div. (35->9), fcn. (880->13), ass. (0->34)
t59 = sin(pkin(10));
t60 = sin(pkin(6));
t71 = t59 * t60;
t62 = cos(pkin(10));
t58 = sin(pkin(11));
t61 = cos(pkin(11));
t65 = sin(qJ(2));
t67 = cos(qJ(2));
t55 = t58 * t65 - t67 * t61;
t63 = cos(pkin(6));
t68 = t55 * t63;
t69 = t58 * t67 + t65 * t61;
t48 = t59 * t68 - t62 * t69;
t64 = sin(qJ(5));
t66 = cos(qJ(5));
t43 = -t48 * t64 + t66 * t71;
t41 = 0.1e1 / t43 ^ 2;
t42 = t48 * t66 + t64 * t71;
t70 = t41 * t42 ^ 2 + 0.1e1;
t54 = t69 * t63;
t49 = -t54 * t59 - t55 * t62;
t53 = t69 * t60;
t52 = t55 * t60;
t51 = 0.1e1 / t53 ^ 2;
t45 = t54 * t62 - t55 * t59;
t44 = -t59 * t69 - t62 * t68;
t40 = atan2(-t45, t53);
t38 = cos(t40);
t37 = sin(t40);
t36 = 0.1e1 / t70;
t35 = -t37 * t45 + t38 * t53;
t34 = 0.1e1 / t35 ^ 2;
t32 = (-t44 / t53 - t52 * t45 * t51) / (t45 ^ 2 * t51 + 0.1e1);
t1 = [0, t32, 0, 0, 0, 0; 0 (t48 / t35 - (-t37 * t44 - t38 * t52 + (-t37 * t53 - t38 * t45) * t32) * t49 * t34) / (t34 * t49 ^ 2 + 0.1e1) 0, 0, 0, 0; 0 (-t66 / t43 - t64 * t42 * t41) * t49 * t36, 0, 0, t70 * t36, 0;];
Ja_rot  = t1;
