% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:18:10
% EndTime: 2019-02-26 21:18:11
% DurationCPUTime: 0.12s
% Computational Cost: add. (694->20), mult. (440->54), div. (106->9), fcn. (646->9), ass. (0->38)
t69 = qJ(3) + qJ(4) + qJ(5);
t67 = sin(t69);
t68 = cos(t69);
t73 = cos(qJ(1));
t77 = t73 * t68;
t62 = atan2(-t77, t67);
t58 = sin(t62);
t59 = cos(t62);
t53 = -t58 * t77 + t59 * t67;
t52 = 0.1e1 / t53 ^ 2;
t71 = sin(qJ(1));
t85 = t52 * t71 ^ 2;
t70 = sin(qJ(6));
t76 = t73 * t70;
t72 = cos(qJ(6));
t79 = t71 * t72;
t61 = t67 * t79 + t76;
t57 = 0.1e1 / t61 ^ 2;
t75 = t73 * t72;
t80 = t71 * t70;
t60 = t67 * t80 - t75;
t84 = t57 * t60;
t83 = t58 * t67;
t66 = t68 ^ 2;
t82 = 0.1e1 / t67 ^ 2 * t66;
t81 = t68 * t71;
t63 = 0.1e1 / (t73 ^ 2 * t82 + 0.1e1);
t78 = t73 * t63;
t74 = t60 ^ 2 * t57 + 0.1e1;
t64 = 0.1e1 / t67;
t56 = 0.1e1 / t61;
t55 = 0.1e1 / t74;
t54 = (0.1e1 + t82) * t78;
t51 = 0.1e1 / t53;
t50 = 0.1e1 / (t66 * t85 + 0.1e1);
t49 = (t56 * t70 - t72 * t84) * t55 * t81;
t48 = (t67 * t51 + (t73 * t83 + t59 * t68 + (-t59 * t77 - t83) * t54) * t68 * t52) * t71 * t50;
t1 = [t64 * t63 * t81, 0, t54, t54, t54, 0; (-t51 * t77 + (-t59 * t64 * t66 * t78 + (-t63 + 0.1e1) * t68 * t58) * t68 * t85) * t50, 0, t48, t48, t48, 0; ((t67 * t76 + t79) * t56 - (t67 * t75 - t80) * t84) * t55, 0, t49, t49, t49, t74 * t55;];
Ja_rot  = t1;
