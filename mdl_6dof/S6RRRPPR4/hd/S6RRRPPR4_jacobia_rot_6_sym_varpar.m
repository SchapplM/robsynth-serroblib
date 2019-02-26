% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:09
% EndTime: 2019-02-26 22:05:09
% DurationCPUTime: 0.10s
% Computational Cost: add. (273->26), mult. (425->68), div. (58->9), fcn. (611->11), ass. (0->41)
t63 = qJ(3) + pkin(10);
t62 = cos(t63);
t69 = sin(qJ(1));
t71 = cos(qJ(2));
t61 = sin(t63);
t72 = cos(qJ(1));
t75 = t72 * t61;
t54 = -t69 * t62 + t71 * t75;
t74 = t72 * t62;
t55 = t69 * t61 + t71 * t74;
t67 = sin(qJ(6));
t70 = cos(qJ(6));
t46 = t54 * t67 + t55 * t70;
t44 = 0.1e1 / t46 ^ 2;
t45 = -t54 * t70 + t55 * t67;
t83 = t44 * t45;
t82 = t45 ^ 2 * t44;
t68 = sin(qJ(2));
t77 = t69 * t68;
t59 = atan2(t77, t71);
t56 = sin(t59);
t57 = cos(t59);
t50 = t56 * t77 + t57 * t71;
t49 = 0.1e1 / t50 ^ 2;
t81 = t49 * t72 ^ 2;
t64 = t68 ^ 2;
t80 = t64 / t71 ^ 2;
t79 = t68 * t72;
t58 = 0.1e1 / (t69 ^ 2 * t80 + 0.1e1);
t78 = t69 * t58;
t76 = t69 * t71;
t73 = 0.1e1 + t82;
t65 = 0.1e1 / t71;
t53 = -t62 * t76 + t75;
t52 = -t61 * t76 - t74;
t51 = (0.1e1 + t80) * t78;
t48 = 0.1e1 / t50;
t47 = 0.1e1 / (t64 * t81 + 0.1e1);
t43 = 0.1e1 / t46;
t42 = 0.1e1 / t73;
t1 = [t65 * t58 * t79, t51, 0, 0, 0, 0; (t48 * t77 + (t57 * t64 * t65 * t78 + (-t58 + 0.1e1) * t68 * t56) * t68 * t81) * t47 (-t71 * t48 + (t56 * t76 - t57 * t68 + (-t56 * t71 + t57 * t77) * t51) * t68 * t49) * t72 * t47, 0, 0, 0, 0; ((-t52 * t70 + t53 * t67) * t43 - (t52 * t67 + t53 * t70) * t83) * t42 ((t61 * t70 - t62 * t67) * t43 - (-t61 * t67 - t62 * t70) * t83) * t42 * t79 (-t46 * t43 - t82) * t42, 0, 0, t73 * t42;];
Ja_rot  = t1;
