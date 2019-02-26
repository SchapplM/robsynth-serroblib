% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRP2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:40:10
% EndTime: 2019-02-26 22:40:10
% DurationCPUTime: 0.09s
% Computational Cost: add. (741->21), mult. (440->54), div. (106->9), fcn. (646->9), ass. (0->38)
t74 = qJ(2) + qJ(3) + qJ(4);
t73 = cos(t74);
t72 = sin(t74);
t76 = sin(qJ(1));
t84 = t76 * t72;
t67 = atan2(-t84, -t73);
t61 = sin(t67);
t62 = cos(t67);
t58 = -t61 * t84 - t62 * t73;
t57 = 0.1e1 / t58 ^ 2;
t78 = cos(qJ(1));
t90 = t57 * t78 ^ 2;
t89 = t61 * t73;
t77 = cos(qJ(5));
t80 = t78 * t77;
t75 = sin(qJ(5));
t83 = t76 * t75;
t66 = t73 * t80 + t83;
t64 = 0.1e1 / t66 ^ 2;
t81 = t78 * t75;
t82 = t76 * t77;
t65 = t73 * t81 - t82;
t88 = t64 * t65;
t69 = t72 ^ 2;
t87 = t69 / t73 ^ 2;
t86 = t72 * t78;
t68 = 0.1e1 / (t76 ^ 2 * t87 + 0.1e1);
t85 = t76 * t68;
t79 = t65 ^ 2 * t64 + 0.1e1;
t70 = 0.1e1 / t73;
t63 = 0.1e1 / t66;
t60 = 0.1e1 / t79;
t59 = (0.1e1 + t87) * t85;
t56 = 0.1e1 / t58;
t55 = 0.1e1 / (t69 * t90 + 0.1e1);
t54 = (-t63 * t75 + t77 * t88) * t60 * t86;
t53 = (t73 * t56 - (-t76 * t89 + t62 * t72 + (-t62 * t84 + t89) * t59) * t72 * t57) * t78 * t55;
t1 = [t70 * t68 * t86, t59, t59, t59, 0, 0; (-t56 * t84 - (-t62 * t69 * t70 * t85 + (t68 - 0.1e1) * t72 * t61) * t72 * t90) * t55, t53, t53, t53, 0, 0; ((-t73 * t83 - t80) * t63 - (-t73 * t82 + t81) * t88) * t60, t54, t54, t54, t79 * t60, 0;];
Ja_rot  = t1;
