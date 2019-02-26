% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:37
% EndTime: 2019-02-26 21:02:37
% DurationCPUTime: 0.11s
% Computational Cost: add. (618->22), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->39)
t70 = pkin(10) + qJ(3) + qJ(4);
t67 = cos(t70);
t66 = sin(t70);
t72 = sin(qJ(1));
t79 = t72 * t66;
t61 = atan2(-t79, -t67);
t59 = sin(t61);
t60 = cos(t61);
t52 = -t59 * t79 - t60 * t67;
t51 = 0.1e1 / t52 ^ 2;
t73 = cos(qJ(1));
t85 = t51 * t73 ^ 2;
t71 = pkin(11) + qJ(6);
t69 = cos(t71);
t75 = t73 * t69;
t68 = sin(t71);
t78 = t72 * t68;
t58 = t67 * t75 + t78;
t56 = 0.1e1 / t58 ^ 2;
t76 = t73 * t68;
t77 = t72 * t69;
t57 = t67 * t76 - t77;
t84 = t56 * t57;
t83 = t59 * t67;
t63 = t66 ^ 2;
t82 = t63 / t67 ^ 2;
t81 = t66 * t73;
t62 = 0.1e1 / (t72 ^ 2 * t82 + 0.1e1);
t80 = t72 * t62;
t74 = t57 ^ 2 * t56 + 0.1e1;
t64 = 0.1e1 / t67;
t55 = 0.1e1 / t58;
t54 = 0.1e1 / t74;
t53 = (0.1e1 + t82) * t80;
t50 = 0.1e1 / t52;
t49 = 0.1e1 / (t63 * t85 + 0.1e1);
t48 = (-t55 * t68 + t69 * t84) * t54 * t81;
t47 = (t67 * t50 - (-t72 * t83 + t60 * t66 + (-t60 * t79 + t83) * t53) * t66 * t51) * t73 * t49;
t1 = [t64 * t62 * t81, 0, t53, t53, 0, 0; (-t50 * t79 - (-t60 * t63 * t64 * t80 + (t62 - 0.1e1) * t66 * t59) * t66 * t85) * t49, 0, t47, t47, 0, 0; ((-t67 * t78 - t75) * t55 - (-t67 * t77 + t76) * t84) * t54, 0, t48, t48, 0, t74 * t54;];
Ja_rot  = t1;
