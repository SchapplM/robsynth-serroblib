% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:28
% EndTime: 2019-02-26 22:03:28
% DurationCPUTime: 0.08s
% Computational Cost: add. (530->21), mult. (305->53), div. (74->9), fcn. (454->9), ass. (0->37)
t67 = qJ(2) + qJ(3) + pkin(10);
t66 = cos(t67);
t65 = sin(t67);
t70 = sin(qJ(1));
t76 = t70 * t65;
t60 = atan2(-t76, -t66);
t54 = sin(t60);
t55 = cos(t60);
t51 = -t54 * t76 - t55 * t66;
t50 = 0.1e1 / t51 ^ 2;
t71 = cos(qJ(1));
t82 = t50 * t71 ^ 2;
t81 = t54 * t66;
t69 = cos(pkin(11));
t72 = t71 * t69;
t68 = sin(pkin(11));
t75 = t70 * t68;
t59 = t66 * t72 + t75;
t57 = 0.1e1 / t59 ^ 2;
t73 = t71 * t68;
t74 = t70 * t69;
t58 = t66 * t73 - t74;
t80 = t57 * t58;
t62 = t65 ^ 2;
t79 = t62 / t66 ^ 2;
t78 = t65 * t71;
t61 = 0.1e1 / (t70 ^ 2 * t79 + 0.1e1);
t77 = t70 * t61;
t63 = 0.1e1 / t66;
t56 = 0.1e1 / t59;
t53 = 0.1e1 / (t58 ^ 2 * t57 + 0.1e1);
t52 = (0.1e1 + t79) * t77;
t49 = 0.1e1 / t51;
t48 = 0.1e1 / (t62 * t82 + 0.1e1);
t47 = (-t56 * t68 + t69 * t80) * t53 * t78;
t46 = (t66 * t49 - (-t70 * t81 + t55 * t65 + (-t55 * t76 + t81) * t52) * t65 * t50) * t71 * t48;
t1 = [t63 * t61 * t78, t52, t52, 0, 0, 0; (-t49 * t76 - (-t55 * t62 * t63 * t77 + (t61 - 0.1e1) * t65 * t54) * t65 * t82) * t48, t46, t46, 0, 0, 0; ((-t66 * t75 - t72) * t56 - (-t66 * t74 + t73) * t80) * t53, t47, t47, 0, 0, 0;];
Ja_rot  = t1;
