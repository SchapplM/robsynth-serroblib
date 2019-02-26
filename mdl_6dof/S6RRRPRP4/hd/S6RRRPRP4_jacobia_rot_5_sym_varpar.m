% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_jacobia_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:11:05
% EndTime: 2019-02-26 22:11:05
% DurationCPUTime: 0.08s
% Computational Cost: add. (324->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
t70 = qJ(2) + qJ(3);
t68 = sin(t70);
t69 = cos(t70);
t72 = sin(qJ(1));
t80 = t72 * t69;
t63 = atan2(-t80, t68);
t61 = sin(t63);
t62 = cos(t63);
t54 = -t61 * t80 + t62 * t68;
t53 = 0.1e1 / t54 ^ 2;
t74 = cos(qJ(1));
t86 = t53 * t74 ^ 2;
t71 = sin(qJ(5));
t77 = t74 * t71;
t73 = cos(qJ(5));
t78 = t72 * t73;
t60 = t68 * t77 + t78;
t58 = 0.1e1 / t60 ^ 2;
t76 = t74 * t73;
t79 = t72 * t71;
t59 = -t68 * t76 + t79;
t85 = t58 * t59;
t84 = t61 * t68;
t67 = t69 ^ 2;
t83 = 0.1e1 / t68 ^ 2 * t67;
t82 = t69 * t74;
t64 = 0.1e1 / (t72 ^ 2 * t83 + 0.1e1);
t81 = t72 * t64;
t75 = t59 ^ 2 * t58 + 0.1e1;
t65 = 0.1e1 / t68;
t57 = 0.1e1 / t60;
t56 = 0.1e1 / t75;
t55 = (0.1e1 + t83) * t81;
t52 = 0.1e1 / t54;
t51 = 0.1e1 / (t67 * t86 + 0.1e1);
t50 = (-t57 * t73 - t71 * t85) * t56 * t82;
t49 = (-t68 * t52 - (t72 * t84 + t62 * t69 + (-t62 * t80 - t84) * t55) * t69 * t53) * t74 * t51;
t1 = [-t65 * t64 * t82, t55, t55, 0, 0, 0; (-t52 * t80 - (t62 * t65 * t67 * t81 + (t64 - 0.1e1) * t69 * t61) * t69 * t86) * t51, t49, t49, 0, 0, 0; ((t68 * t78 + t77) * t57 - (-t68 * t79 + t76) * t85) * t56, t50, t50, 0, t75 * t56, 0;];
Ja_rot  = t1;
