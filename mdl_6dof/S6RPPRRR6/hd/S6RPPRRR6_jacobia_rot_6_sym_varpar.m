% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:37:35
% EndTime: 2019-02-26 20:37:35
% DurationCPUTime: 0.10s
% Computational Cost: add. (135->19), mult. (251->54), div. (57->9), fcn. (367->9), ass. (0->37)
t54 = sin(qJ(4));
t55 = sin(qJ(1));
t56 = cos(qJ(4));
t62 = t55 * t56;
t47 = atan2(t62, t54);
t44 = sin(t47);
t45 = cos(t47);
t38 = t44 * t62 + t45 * t54;
t37 = 0.1e1 / t38 ^ 2;
t57 = cos(qJ(1));
t69 = t37 * t57 ^ 2;
t53 = qJ(5) + qJ(6);
t49 = cos(t53);
t59 = t57 * t49;
t48 = sin(t53);
t64 = t55 * t48;
t43 = t54 * t59 - t64;
t41 = 0.1e1 / t43 ^ 2;
t60 = t57 * t48;
t63 = t55 * t49;
t42 = t54 * t60 + t63;
t68 = t41 * t42;
t67 = t44 * t54;
t52 = t56 ^ 2;
t66 = 0.1e1 / t54 ^ 2 * t52;
t46 = 0.1e1 / (t55 ^ 2 * t66 + 0.1e1);
t65 = t55 * t46;
t61 = t56 * t57;
t58 = t42 ^ 2 * t41 + 0.1e1;
t50 = 0.1e1 / t54;
t40 = 0.1e1 / t43;
t39 = (-0.1e1 - t66) * t65;
t36 = 0.1e1 / t38;
t35 = 0.1e1 / t58;
t34 = 0.1e1 / (t52 * t69 + 0.1e1);
t33 = t58 * t35;
t1 = [t50 * t46 * t61, 0, 0, t39, 0, 0; (t36 * t62 + (t45 * t50 * t52 * t65 + (-t46 + 0.1e1) * t56 * t44) * t56 * t69) * t34, 0, 0 (t54 * t36 + (-t55 * t67 + t45 * t56 + (t45 * t62 - t67) * t39) * t56 * t37) * t57 * t34, 0, 0; ((-t54 * t64 + t59) * t40 - (-t54 * t63 - t60) * t68) * t35, 0, 0 (t40 * t48 - t49 * t68) * t35 * t61, t33, t33;];
Ja_rot  = t1;
