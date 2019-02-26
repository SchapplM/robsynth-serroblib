% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR5
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
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:48
% EndTime: 2019-02-26 21:30:48
% DurationCPUTime: 0.08s
% Computational Cost: add. (149->22), mult. (451->62), div. (72->11), fcn. (673->11), ass. (0->42)
t50 = cos(pkin(6));
t55 = cos(qJ(2));
t56 = cos(qJ(1));
t58 = t56 * t55;
t52 = sin(qJ(2));
t53 = sin(qJ(1));
t61 = t53 * t52;
t44 = -t50 * t61 + t58;
t51 = sin(qJ(5));
t54 = cos(qJ(5));
t49 = sin(pkin(6));
t64 = t49 * t53;
t35 = t44 * t54 - t51 * t64;
t33 = 0.1e1 / t35 ^ 2;
t34 = t44 * t51 + t54 * t64;
t68 = t33 * t34;
t40 = -t50 * t58 + t61;
t63 = t49 * t55;
t39 = atan2(t40, t63);
t37 = cos(t39);
t67 = t37 * t40;
t36 = sin(t39);
t30 = t36 * t40 + t37 * t63;
t29 = 0.1e1 / t30 ^ 2;
t59 = t56 * t52;
t60 = t53 * t55;
t42 = t50 * t60 + t59;
t66 = t42 ^ 2 * t29;
t46 = 0.1e1 / t49;
t47 = 0.1e1 / t55;
t65 = t46 * t47;
t62 = t49 * t56;
t57 = t34 ^ 2 * t33 + 0.1e1;
t48 = 0.1e1 / t55 ^ 2;
t41 = t50 * t59 + t60;
t38 = 0.1e1 / (0.1e1 + t40 ^ 2 / t49 ^ 2 * t48);
t32 = 0.1e1 / t35;
t31 = 0.1e1 / t57;
t28 = 0.1e1 / t30;
t27 = 0.1e1 / (0.1e1 + t66);
t26 = (t40 * t48 * t52 + t41 * t47) * t46 * t38;
t1 = [t42 * t38 * t65, t26, 0, 0, 0, 0; (t40 * t28 + (t36 + (t65 * t67 - t36) * t38) * t66) * t27 (-t44 * t28 + (-t37 * t49 * t52 + t36 * t41 + (-t36 * t63 + t67) * t26) * t42 * t29) * t27, 0, 0, 0, 0; ((-t41 * t51 + t54 * t62) * t32 - (-t41 * t54 - t51 * t62) * t68) * t31 (-t51 * t32 + t54 * t68) * t42 * t31, 0, 0, t57 * t31, 0;];
Ja_rot  = t1;
