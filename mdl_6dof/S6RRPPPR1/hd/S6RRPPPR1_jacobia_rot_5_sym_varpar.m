% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPPR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:58
% EndTime: 2019-02-26 21:21:59
% DurationCPUTime: 0.10s
% Computational Cost: add. (227->21), mult. (329->51), div. (61->11), fcn. (494->9), ass. (0->36)
t52 = qJ(2) + pkin(9);
t49 = sin(t52);
t69 = t49 ^ 2;
t50 = cos(t52);
t55 = cos(pkin(10));
t57 = cos(qJ(1));
t59 = t57 * t55;
t54 = sin(pkin(10));
t56 = sin(qJ(1));
t62 = t56 * t54;
t41 = t50 * t62 + t59;
t63 = t49 * t54;
t37 = atan2(-t41, t63);
t34 = sin(t37);
t35 = cos(t37);
t33 = -t34 * t41 + t35 * t63;
t32 = 0.1e1 / t33 ^ 2;
t60 = t57 * t54;
t61 = t56 * t55;
t43 = t50 * t60 - t61;
t68 = t32 * t43;
t66 = t35 * t41;
t65 = t43 ^ 2 * t32;
t51 = 0.1e1 / t54;
t64 = 0.1e1 / t49 * t51;
t44 = t50 * t59 + t62;
t40 = 0.1e1 / t44 ^ 2;
t58 = t57 ^ 2 * t69 * t40;
t48 = 0.1e1 / t69;
t39 = 0.1e1 / t44;
t38 = 0.1e1 / (0.1e1 + t58);
t36 = 0.1e1 / (0.1e1 + t41 ^ 2 * t48 / t54 ^ 2);
t31 = 0.1e1 / t33;
t30 = (t41 * t48 * t50 * t51 + t56) * t36;
t29 = 0.1e1 / (0.1e1 + t65);
t1 = [-t43 * t36 * t64, t30, 0, 0, 0, 0; (-t41 * t31 - (-t34 + (t64 * t66 + t34) * t36) * t65) * t29 (t30 * t66 * t68 + (-t57 * t49 * t31 - (t35 * t50 + (-t30 + t56) * t34 * t49) * t68) * t54) * t29, 0, 0, 0, 0; (t56 * t39 + (-t50 * t61 + t60) * t57 * t40) * t49 * t38 (-t39 * t50 * t57 - t55 * t58) * t38, 0, 0, 0, 0;];
Ja_rot  = t1;
