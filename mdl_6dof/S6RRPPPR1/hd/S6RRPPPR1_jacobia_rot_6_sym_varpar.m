% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_rot = S6RRPPPR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:49
% EndTime: 2019-02-26 21:21:49
% DurationCPUTime: 0.12s
% Computational Cost: add. (265->25), mult. (347->66), div. (52->9), fcn. (503->11), ass. (0->42)
t58 = qJ(2) + pkin(9);
t57 = cos(t58);
t56 = sin(t58);
t62 = sin(qJ(1));
t70 = t62 * t56;
t52 = atan2(t70, t57);
t49 = sin(t52);
t50 = cos(t52);
t39 = t49 * t70 + t50 * t57;
t38 = 0.1e1 / t39 ^ 2;
t64 = cos(qJ(1));
t76 = t38 * t64 ^ 2;
t59 = sin(pkin(10));
t67 = t64 * t59;
t60 = cos(pkin(10));
t68 = t62 * t60;
t47 = t57 * t67 - t68;
t66 = t64 * t60;
t69 = t62 * t59;
t48 = t57 * t66 + t69;
t61 = sin(qJ(6));
t63 = cos(qJ(6));
t44 = t47 * t61 + t48 * t63;
t42 = 0.1e1 / t44 ^ 2;
t43 = -t47 * t63 + t48 * t61;
t75 = t42 * t43;
t74 = t49 * t57;
t53 = t56 ^ 2;
t73 = t53 / t57 ^ 2;
t72 = t56 * t64;
t51 = 0.1e1 / (t62 ^ 2 * t73 + 0.1e1);
t71 = t62 * t51;
t65 = t42 * t43 ^ 2 + 0.1e1;
t54 = 0.1e1 / t57;
t46 = -t57 * t68 + t67;
t45 = -t57 * t69 - t66;
t41 = 0.1e1 / t44;
t40 = (0.1e1 + t73) * t71;
t37 = 0.1e1 / t39;
t36 = 0.1e1 / (t53 * t76 + 0.1e1);
t35 = 0.1e1 / t65;
t1 = [t54 * t51 * t72, t40, 0, 0, 0, 0; (t37 * t70 + (t50 * t53 * t54 * t71 + (-t51 + 0.1e1) * t56 * t49) * t56 * t76) * t36 (-t57 * t37 + (t62 * t74 - t50 * t56 + (t50 * t70 - t74) * t40) * t56 * t38) * t64 * t36, 0, 0, 0, 0; ((-t45 * t63 + t46 * t61) * t41 - (t45 * t61 + t46 * t63) * t75) * t35 ((t59 * t63 - t60 * t61) * t41 - (-t59 * t61 - t60 * t63) * t75) * t35 * t72, 0, 0, 0, t65 * t35;];
Ja_rot  = t1;
