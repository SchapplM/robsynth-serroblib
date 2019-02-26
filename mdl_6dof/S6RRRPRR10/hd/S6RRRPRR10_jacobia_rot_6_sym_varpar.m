% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR10_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:23
% EndTime: 2019-02-26 22:21:23
% DurationCPUTime: 0.10s
% Computational Cost: add. (247->26), mult. (487->68), div. (63->9), fcn. (695->11), ass. (0->42)
t77 = sin(qJ(1));
t78 = cos(qJ(3));
t79 = cos(qJ(2));
t75 = sin(qJ(3));
t80 = cos(qJ(1));
t83 = t80 * t75;
t62 = -t77 * t78 + t79 * t83;
t82 = t80 * t78;
t63 = t77 * t75 + t79 * t82;
t74 = qJ(5) + qJ(6);
t69 = sin(t74);
t70 = cos(t74);
t54 = t62 * t69 + t63 * t70;
t52 = 0.1e1 / t54 ^ 2;
t53 = -t62 * t70 + t63 * t69;
t91 = t52 * t53;
t90 = t53 ^ 2 * t52;
t76 = sin(qJ(2));
t85 = t77 * t76;
t67 = atan2(t85, t79);
t64 = sin(t67);
t65 = cos(t67);
t58 = t64 * t85 + t65 * t79;
t57 = 0.1e1 / t58 ^ 2;
t89 = t57 * t80 ^ 2;
t71 = t76 ^ 2;
t88 = t71 / t79 ^ 2;
t87 = t76 * t80;
t66 = 0.1e1 / (t77 ^ 2 * t88 + 0.1e1);
t86 = t77 * t66;
t84 = t77 * t79;
t81 = 0.1e1 + t90;
t72 = 0.1e1 / t79;
t61 = -t78 * t84 + t83;
t60 = -t75 * t84 - t82;
t59 = (0.1e1 + t88) * t86;
t56 = 0.1e1 / t58;
t55 = 0.1e1 / (t71 * t89 + 0.1e1);
t51 = 0.1e1 / t54;
t50 = 0.1e1 / t81;
t49 = t81 * t50;
t1 = [t72 * t66 * t87, t59, 0, 0, 0, 0; (t56 * t85 + (t65 * t71 * t72 * t86 + (-t66 + 0.1e1) * t76 * t64) * t76 * t89) * t55 (-t79 * t56 + (t64 * t84 - t65 * t76 + (-t64 * t79 + t65 * t85) * t59) * t76 * t57) * t80 * t55, 0, 0, 0, 0; ((-t60 * t70 + t61 * t69) * t51 - (t60 * t69 + t61 * t70) * t91) * t50 ((-t69 * t78 + t70 * t75) * t51 - (-t69 * t75 - t70 * t78) * t91) * t50 * t87 (-t54 * t51 - t90) * t50, 0, t49, t49;];
Ja_rot  = t1;
