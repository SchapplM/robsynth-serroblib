% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:18
% EndTime: 2019-02-26 22:03:18
% DurationCPUTime: 0.08s
% Computational Cost: add. (618->22), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->39)
t75 = qJ(2) + qJ(3) + pkin(10);
t72 = cos(t75);
t71 = sin(t75);
t77 = sin(qJ(1));
t84 = t77 * t71;
t66 = atan2(-t84, -t72);
t64 = sin(t66);
t65 = cos(t66);
t57 = -t64 * t84 - t65 * t72;
t56 = 0.1e1 / t57 ^ 2;
t78 = cos(qJ(1));
t90 = t56 * t78 ^ 2;
t76 = pkin(11) + qJ(6);
t74 = cos(t76);
t80 = t78 * t74;
t73 = sin(t76);
t83 = t77 * t73;
t63 = t72 * t80 + t83;
t61 = 0.1e1 / t63 ^ 2;
t81 = t78 * t73;
t82 = t77 * t74;
t62 = t72 * t81 - t82;
t89 = t61 * t62;
t88 = t64 * t72;
t68 = t71 ^ 2;
t87 = t68 / t72 ^ 2;
t86 = t71 * t78;
t67 = 0.1e1 / (t77 ^ 2 * t87 + 0.1e1);
t85 = t77 * t67;
t79 = t62 ^ 2 * t61 + 0.1e1;
t69 = 0.1e1 / t72;
t60 = 0.1e1 / t63;
t59 = 0.1e1 / t79;
t58 = (0.1e1 + t87) * t85;
t55 = 0.1e1 / t57;
t54 = 0.1e1 / (t68 * t90 + 0.1e1);
t53 = (-t60 * t73 + t74 * t89) * t59 * t86;
t52 = (t72 * t55 - (-t77 * t88 + t65 * t71 + (-t65 * t84 + t88) * t58) * t71 * t56) * t78 * t54;
t1 = [t69 * t67 * t86, t58, t58, 0, 0, 0; (-t55 * t84 - (-t65 * t68 * t69 * t85 + (t67 - 0.1e1) * t71 * t64) * t71 * t90) * t54, t52, t52, 0, 0, 0; ((-t72 * t83 - t80) * t60 - (-t72 * t82 + t81) * t89) * t59, t53, t53, 0, 0, t79 * t59;];
Ja_rot  = t1;
