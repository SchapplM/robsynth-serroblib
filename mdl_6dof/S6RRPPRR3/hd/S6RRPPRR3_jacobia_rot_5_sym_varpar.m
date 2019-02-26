% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_rot = S6RRPPRR3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:33
% EndTime: 2019-02-26 21:29:33
% DurationCPUTime: 0.18s
% Computational Cost: add. (422->27), mult. (1051->69), div. (55->9), fcn. (1490->13), ass. (0->46)
t82 = sin(qJ(1));
t84 = cos(qJ(1));
t78 = sin(pkin(11));
t80 = cos(pkin(11));
t83 = cos(qJ(2));
t91 = cos(pkin(6));
t88 = t83 * t91;
t81 = sin(qJ(2));
t89 = t81 * t91;
t85 = -t78 * t89 + t80 * t88;
t86 = t83 * t78 + t81 * t80;
t61 = -t82 * t86 + t84 * t85;
t72 = t81 * t78 - t83 * t80;
t79 = sin(pkin(6));
t69 = t72 * t79;
t55 = atan2(t61, t69);
t53 = cos(t55);
t96 = t53 * t61;
t71 = t78 * t88 + t80 * t89;
t65 = -t82 * t71 - t84 * t72;
t77 = pkin(12) + qJ(5);
t75 = sin(t77);
t76 = cos(t77);
t93 = t79 * t82;
t59 = t65 * t76 + t75 * t93;
t57 = 0.1e1 / t59 ^ 2;
t58 = t65 * t75 - t76 * t93;
t95 = t57 * t58;
t52 = sin(t55);
t50 = t52 * t61 + t53 * t69;
t49 = 0.1e1 / t50 ^ 2;
t63 = -t82 * t85 - t84 * t86;
t94 = t63 ^ 2 * t49;
t92 = t79 * t84;
t90 = t58 ^ 2 * t57 + 0.1e1;
t87 = -t84 * t71 + t82 * t72;
t70 = t86 * t79;
t68 = 0.1e1 / t69 ^ 2;
t67 = 0.1e1 / t69;
t56 = 0.1e1 / t59;
t54 = 0.1e1 / (t61 ^ 2 * t68 + 0.1e1);
t51 = 0.1e1 / t90;
t48 = 0.1e1 / t50;
t47 = 0.1e1 / (0.1e1 + t94);
t46 = (-t61 * t68 * t70 + t67 * t87) * t54;
t1 = [t63 * t67 * t54, t46, 0, 0, 0, 0; (t61 * t48 + (t52 + (t67 * t96 - t52) * t54) * t94) * t47 (t65 * t48 + (t52 * t87 + t53 * t70 + (-t52 * t69 + t96) * t46) * t63 * t49) * t47, 0, 0, 0, 0; ((t75 * t87 - t76 * t92) * t56 - (t75 * t92 + t76 * t87) * t95) * t51 (t75 * t56 - t76 * t95) * t63 * t51, 0, 0, t90 * t51, 0;];
Ja_rot  = t1;
