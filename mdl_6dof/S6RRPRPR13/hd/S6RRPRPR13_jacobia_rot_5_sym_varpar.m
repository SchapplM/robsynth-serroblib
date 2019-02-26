% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR13_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:50
% EndTime: 2019-02-26 21:44:50
% DurationCPUTime: 0.16s
% Computational Cost: add. (428->36), mult. (1217->90), div. (80->9), fcn. (1746->13), ass. (0->53)
t69 = cos(pkin(6));
t74 = cos(qJ(2));
t75 = cos(qJ(1));
t79 = t75 * t74;
t71 = sin(qJ(2));
t72 = sin(qJ(1));
t82 = t72 * t71;
t62 = -t69 * t79 + t82;
t70 = sin(qJ(4));
t73 = cos(qJ(4));
t67 = sin(pkin(6));
t83 = t67 * t75;
t54 = t62 * t73 + t70 * t83;
t84 = t67 * t74;
t60 = t69 * t70 + t73 * t84;
t51 = atan2(t54, t60);
t48 = sin(t51);
t49 = cos(t51);
t42 = t48 * t54 + t49 * t60;
t41 = 0.1e1 / t42 ^ 2;
t80 = t75 * t71;
t81 = t72 * t74;
t76 = t69 * t81 + t80;
t85 = t67 * t72;
t52 = t70 * t85 - t73 * t76;
t92 = t41 * t52;
t53 = t70 * t76 + t73 * t85;
t64 = -t69 * t82 + t79;
t66 = sin(pkin(11));
t68 = cos(pkin(11));
t47 = t53 * t68 + t64 * t66;
t45 = 0.1e1 / t47 ^ 2;
t46 = t53 * t66 - t64 * t68;
t91 = t45 * t46;
t90 = t49 * t54;
t89 = t52 ^ 2 * t41;
t59 = 0.1e1 / t60 ^ 2;
t88 = t54 * t59;
t87 = t64 * t70;
t86 = t67 * t71;
t78 = -t48 * t60 + t90;
t77 = -t62 * t70 + t73 * t83;
t63 = t69 * t80 + t81;
t61 = t69 * t73 - t70 * t84;
t58 = 0.1e1 / t60;
t50 = 0.1e1 / (t54 ^ 2 * t59 + 0.1e1);
t44 = 0.1e1 / t47;
t43 = 0.1e1 / (t45 * t46 ^ 2 + 0.1e1);
t40 = 0.1e1 / t42;
t39 = 0.1e1 / (0.1e1 + t89);
t38 = (t58 * t63 + t86 * t88) * t73 * t50;
t37 = (t58 * t77 - t61 * t88) * t50;
t1 = [-t52 * t58 * t50, t38, 0, t37, 0, 0; (t54 * t40 - (-t48 + (-t58 * t90 + t48) * t50) * t89) * t39 (-t64 * t73 * t40 - ((t48 * t63 - t49 * t86) * t73 + t78 * t38) * t92) * t39, 0 (t53 * t40 - (t37 * t78 + t48 * t77 + t49 * t61) * t92) * t39, 0, 0; ((t63 * t68 + t66 * t77) * t44 - (-t63 * t66 + t68 * t77) * t91) * t43 ((t66 * t87 + t68 * t76) * t44 - (-t66 * t76 + t68 * t87) * t91) * t43, 0 (-t44 * t66 + t68 * t91) * t52 * t43, 0, 0;];
Ja_rot  = t1;
