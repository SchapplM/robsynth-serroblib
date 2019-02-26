% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR8_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:16
% EndTime: 2019-02-26 21:58:16
% DurationCPUTime: 0.07s
% Computational Cost: add. (357->21), mult. (278->55), div. (62->9), fcn. (402->9), ass. (0->35)
t66 = cos(qJ(2));
t64 = sin(qJ(2));
t65 = sin(qJ(1));
t72 = t65 * t64;
t56 = atan2(-t72, -t66);
t54 = sin(t56);
t55 = cos(t56);
t48 = -t54 * t72 - t55 * t66;
t47 = 0.1e1 / t48 ^ 2;
t67 = cos(qJ(1));
t77 = t47 * t67 ^ 2;
t60 = pkin(11) + qJ(4) + qJ(5) + qJ(6);
t58 = sin(t60);
t59 = cos(t60);
t69 = t67 * t59;
t53 = t65 * t58 + t66 * t69;
t51 = 0.1e1 / t53 ^ 2;
t70 = t67 * t58;
t52 = -t65 * t59 + t66 * t70;
t76 = t51 * t52;
t61 = t64 ^ 2;
t75 = t61 / t66 ^ 2;
t74 = t64 * t67;
t57 = 0.1e1 / (t65 ^ 2 * t75 + 0.1e1);
t73 = t65 * t57;
t71 = t65 * t66;
t68 = t52 ^ 2 * t51 + 0.1e1;
t62 = 0.1e1 / t66;
t50 = 0.1e1 / t53;
t49 = (0.1e1 + t75) * t73;
t46 = 0.1e1 / t48;
t45 = 0.1e1 / (t61 * t77 + 0.1e1);
t44 = 0.1e1 / t68;
t43 = t68 * t44;
t1 = [t62 * t57 * t74, t49, 0, 0, 0, 0; (-t46 * t72 - (-t55 * t61 * t62 * t73 + (t57 - 0.1e1) * t64 * t54) * t64 * t77) * t45 (t66 * t46 - (-t54 * t71 + t55 * t64 + (t54 * t66 - t55 * t72) * t49) * t64 * t47) * t67 * t45, 0, 0, 0, 0; ((-t58 * t71 - t69) * t50 - (-t59 * t71 + t70) * t76) * t44 (-t50 * t58 + t59 * t76) * t44 * t74, 0, t43, t43, t43;];
Ja_rot  = t1;
