% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP12_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:52:13
% EndTime: 2019-02-26 21:52:13
% DurationCPUTime: 0.14s
% Computational Cost: add. (610->27), mult. (689->67), div. (139->11), fcn. (1046->9), ass. (0->43)
t70 = cos(qJ(2));
t84 = t70 ^ 2;
t68 = sin(qJ(2));
t67 = qJ(4) + qJ(5);
t61 = sin(t67);
t71 = cos(qJ(1));
t74 = t71 * t61;
t62 = cos(t67);
t69 = sin(qJ(1));
t77 = t69 * t62;
t56 = t68 * t77 + t74;
t76 = t70 * t62;
t50 = atan2(t56, t76);
t47 = sin(t50);
t48 = cos(t50);
t45 = t47 * t56 + t48 * t76;
t44 = 0.1e1 / t45 ^ 2;
t73 = t71 * t62;
t78 = t69 * t61;
t54 = -t68 * t73 + t78;
t83 = t44 * t54;
t81 = t48 * t56;
t80 = t54 ^ 2 * t44;
t59 = 0.1e1 / t62;
t64 = 0.1e1 / t70;
t79 = t59 * t64;
t75 = t70 * t71;
t55 = t68 * t74 + t77;
t53 = 0.1e1 / t55 ^ 2;
t72 = t71 ^ 2 * t84 * t53;
t65 = 0.1e1 / t84;
t60 = 0.1e1 / t62 ^ 2;
t57 = -t68 * t78 + t73;
t52 = 0.1e1 / t55;
t51 = 0.1e1 / (0.1e1 + t72);
t49 = 0.1e1 / (t56 ^ 2 * t65 * t60 + 0.1e1);
t46 = t54 * t53 * t51 * t75;
t43 = 0.1e1 / t45;
t42 = (t56 * t59 * t65 * t68 + t69) * t49;
t41 = 0.1e1 / (0.1e1 + t80);
t40 = (t56 * t60 * t61 + t57 * t59) * t64 * t49;
t39 = (t55 * t43 - (-t48 * t70 * t61 + t47 * t57 + (-t47 * t76 + t81) * t40) * t83) * t41;
t1 = [-t54 * t49 * t79, t42, 0, t40, t40, 0; (t56 * t43 - (-t47 + (-t79 * t81 + t47) * t49) * t80) * t41 (-t42 * t81 * t83 + (-t43 * t75 - (-t48 * t68 + (-t42 + t69) * t70 * t47) * t83) * t62) * t41, 0, t39, t39, 0; (t53 * t57 * t71 + t52 * t69) * t70 * t51 (t52 * t68 * t71 + t61 * t72) * t51, 0, -t46, -t46, 0;];
Ja_rot  = t1;
