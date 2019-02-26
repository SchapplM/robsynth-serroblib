% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR9_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:33:12
% EndTime: 2019-02-26 21:33:12
% DurationCPUTime: 0.08s
% Computational Cost: add. (149->22), mult. (451->62), div. (72->11), fcn. (673->11), ass. (0->42)
t52 = cos(pkin(6));
t57 = cos(qJ(2));
t58 = cos(qJ(1));
t60 = t58 * t57;
t54 = sin(qJ(2));
t55 = sin(qJ(1));
t63 = t55 * t54;
t46 = -t52 * t63 + t60;
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t51 = sin(pkin(6));
t66 = t51 * t55;
t37 = t46 * t53 + t56 * t66;
t35 = 0.1e1 / t37 ^ 2;
t36 = -t46 * t56 + t53 * t66;
t70 = t35 * t36;
t42 = -t52 * t60 + t63;
t65 = t51 * t57;
t41 = atan2(t42, t65);
t39 = cos(t41);
t69 = t39 * t42;
t38 = sin(t41);
t32 = t38 * t42 + t39 * t65;
t31 = 0.1e1 / t32 ^ 2;
t61 = t58 * t54;
t62 = t55 * t57;
t44 = t52 * t62 + t61;
t68 = t44 ^ 2 * t31;
t48 = 0.1e1 / t51;
t49 = 0.1e1 / t57;
t67 = t48 * t49;
t64 = t51 * t58;
t59 = t36 ^ 2 * t35 + 0.1e1;
t50 = 0.1e1 / t57 ^ 2;
t43 = t52 * t61 + t62;
t40 = 0.1e1 / (0.1e1 + t42 ^ 2 / t51 ^ 2 * t50);
t34 = 0.1e1 / t37;
t33 = 0.1e1 / t59;
t30 = 0.1e1 / t32;
t29 = 0.1e1 / (0.1e1 + t68);
t28 = (t42 * t50 * t54 + t43 * t49) * t48 * t40;
t1 = [t44 * t40 * t67, t28, 0, 0, 0, 0; (t42 * t30 + (t38 + (t67 * t69 - t38) * t40) * t68) * t29 (-t46 * t30 + (-t39 * t51 * t54 + t38 * t43 + (-t38 * t65 + t69) * t28) * t44 * t31) * t29, 0, 0, 0, 0; ((t43 * t56 + t53 * t64) * t34 - (-t43 * t53 + t56 * t64) * t70) * t33 (t56 * t34 + t53 * t70) * t44 * t33, 0, 0, t59 * t33, 0;];
Ja_rot  = t1;
