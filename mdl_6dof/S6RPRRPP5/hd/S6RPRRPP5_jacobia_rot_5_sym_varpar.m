% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:58:48
% EndTime: 2019-02-26 20:58:48
% DurationCPUTime: 0.10s
% Computational Cost: add. (353->27), mult. (509->69), div. (100->11), fcn. (770->9), ass. (0->41)
t53 = pkin(9) + qJ(3);
t51 = sin(t53);
t73 = t51 ^ 2;
t52 = cos(t53);
t59 = cos(qJ(4));
t60 = cos(qJ(1));
t62 = t60 * t59;
t57 = sin(qJ(4));
t58 = sin(qJ(1));
t65 = t58 * t57;
t41 = t52 * t65 + t62;
t67 = t51 * t57;
t38 = atan2(-t41, t67);
t34 = sin(t38);
t35 = cos(t38);
t33 = -t34 * t41 + t35 * t67;
t32 = 0.1e1 / t33 ^ 2;
t63 = t60 * t57;
t64 = t58 * t59;
t44 = t52 * t63 - t64;
t72 = t32 * t44;
t70 = t35 * t41;
t69 = t44 ^ 2 * t32;
t49 = 0.1e1 / t51;
t54 = 0.1e1 / t57;
t68 = t49 * t54;
t66 = t51 * t60;
t45 = t52 * t62 + t65;
t40 = 0.1e1 / t45 ^ 2;
t61 = t60 ^ 2 * t73 * t40;
t55 = 0.1e1 / t57 ^ 2;
t50 = 0.1e1 / t73;
t43 = t52 * t64 - t63;
t39 = 0.1e1 / t45;
t37 = 0.1e1 / (t41 ^ 2 * t50 * t55 + 0.1e1);
t36 = 0.1e1 / (0.1e1 + t61);
t31 = 0.1e1 / t33;
t30 = (t41 * t50 * t52 * t54 + t58) * t37;
t29 = 0.1e1 / (0.1e1 + t69);
t28 = (t41 * t55 * t59 - t43 * t54) * t49 * t37;
t1 = [-t44 * t37 * t68, 0, t30, t28, 0, 0; (-t41 * t31 - (-t34 + (t68 * t70 + t34) * t37) * t69) * t29, 0 (t30 * t70 * t72 + (-t31 * t66 - (t35 * t52 + (-t30 + t58) * t51 * t34) * t72) * t57) * t29 (t45 * t31 - (t35 * t51 * t59 - t34 * t43 + (-t34 * t67 - t70) * t28) * t72) * t29, 0, 0; (-t40 * t43 * t60 + t39 * t58) * t51 * t36, 0 (-t39 * t52 * t60 - t59 * t61) * t36, -t44 * t40 * t36 * t66, 0, 0;];
Ja_rot  = t1;
