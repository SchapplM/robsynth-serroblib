% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:42
% EndTime: 2019-02-26 20:11:43
% DurationCPUTime: 0.12s
% Computational Cost: add. (788->28), mult. (1092->68), div. (84->9), fcn. (1576->11), ass. (0->46)
t78 = sin(pkin(11));
t80 = cos(pkin(11));
t83 = cos(qJ(2));
t81 = cos(pkin(6));
t82 = sin(qJ(2));
t86 = t81 * t82;
t72 = t78 * t83 + t80 * t86;
t77 = qJ(3) + qJ(4);
t75 = sin(t77);
t76 = cos(t77);
t79 = sin(pkin(6));
t89 = t79 * t80;
t60 = t72 * t75 + t76 * t89;
t88 = t79 * t82;
t67 = t75 * t88 - t81 * t76;
t58 = atan2(-t60, t67);
t55 = sin(t58);
t56 = cos(t58);
t53 = -t55 * t60 + t56 * t67;
t52 = 0.1e1 / t53 ^ 2;
t74 = -t78 * t86 + t80 * t83;
t90 = t78 * t79;
t63 = t74 * t75 - t76 * t90;
t92 = t52 * t63;
t66 = 0.1e1 / t67 ^ 2;
t91 = t60 * t66;
t87 = t79 * t83;
t85 = t81 * t83;
t84 = -t55 * t67 - t56 * t60;
t73 = t78 * t85 + t80 * t82;
t71 = -t78 * t82 + t80 * t85;
t70 = 0.1e1 / t73 ^ 2;
t69 = 0.1e1 / t73;
t68 = t81 * t75 + t76 * t88;
t65 = 0.1e1 / t67;
t64 = t74 * t76 + t75 * t90;
t62 = t72 * t76 - t75 * t89;
t59 = 0.1e1 / (t64 ^ 2 * t70 + 0.1e1);
t57 = 0.1e1 / (t60 ^ 2 * t66 + 0.1e1);
t54 = t63 * t69 * t59;
t51 = 0.1e1 / t53;
t50 = 0.1e1 / (t63 ^ 2 * t52 + 0.1e1);
t49 = (-t65 * t71 + t87 * t91) * t75 * t57;
t48 = (-t62 * t65 + t68 * t91) * t57;
t47 = (t64 * t51 - (t84 * t48 - t55 * t62 + t56 * t68) * t92) * t50;
t1 = [0, t49, t48, t48, 0, 0; 0 (-t73 * t75 * t51 - ((-t55 * t71 + t56 * t87) * t75 + t84 * t49) * t92) * t50, t47, t47, 0, 0; 0 (-t64 * t70 * t74 - t69 * t73 * t76) * t59, -t54, -t54, 0, 0;];
Ja_rot  = t1;
