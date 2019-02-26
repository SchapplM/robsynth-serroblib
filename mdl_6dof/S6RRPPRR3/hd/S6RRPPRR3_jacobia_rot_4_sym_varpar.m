% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR3
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR3_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobia_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:33
% EndTime: 2019-02-26 21:29:33
% DurationCPUTime: 0.15s
% Computational Cost: add. (343->26), mult. (968->68), div. (50->9), fcn. (1378->13), ass. (0->44)
t73 = sin(qJ(1));
t75 = cos(qJ(1));
t68 = sin(pkin(11));
t71 = cos(pkin(11));
t74 = cos(qJ(2));
t81 = cos(pkin(6));
t79 = t74 * t81;
t72 = sin(qJ(2));
t80 = t72 * t81;
t76 = -t68 * t80 + t71 * t79;
t77 = t74 * t68 + t72 * t71;
t53 = -t73 * t77 + t75 * t76;
t64 = t72 * t68 - t74 * t71;
t69 = sin(pkin(6));
t61 = t64 * t69;
t47 = atan2(t53, t61);
t45 = cos(t47);
t86 = t45 * t53;
t63 = t68 * t79 + t71 * t80;
t57 = -t73 * t63 - t75 * t64;
t67 = sin(pkin(12));
t70 = cos(pkin(12));
t83 = t69 * t73;
t51 = t57 * t70 + t67 * t83;
t49 = 0.1e1 / t51 ^ 2;
t50 = t57 * t67 - t70 * t83;
t85 = t49 * t50;
t44 = sin(t47);
t42 = t44 * t53 + t45 * t61;
t41 = 0.1e1 / t42 ^ 2;
t55 = -t73 * t76 - t75 * t77;
t84 = t55 ^ 2 * t41;
t82 = t69 * t75;
t78 = -t75 * t63 + t73 * t64;
t62 = t77 * t69;
t60 = 0.1e1 / t61 ^ 2;
t59 = 0.1e1 / t61;
t48 = 0.1e1 / t51;
t46 = 0.1e1 / (t53 ^ 2 * t60 + 0.1e1);
t43 = 0.1e1 / (t50 ^ 2 * t49 + 0.1e1);
t40 = 0.1e1 / t42;
t39 = 0.1e1 / (0.1e1 + t84);
t38 = (-t53 * t60 * t62 + t59 * t78) * t46;
t1 = [t55 * t59 * t46, t38, 0, 0, 0, 0; (t53 * t40 + (t44 + (t59 * t86 - t44) * t46) * t84) * t39 (t57 * t40 + (t44 * t78 + t45 * t62 + (-t44 * t61 + t86) * t38) * t55 * t41) * t39, 0, 0, 0, 0; ((t67 * t78 - t70 * t82) * t48 - (t67 * t82 + t70 * t78) * t85) * t43 (t67 * t48 - t70 * t85) * t55 * t43, 0, 0, 0, 0;];
Ja_rot  = t1;
