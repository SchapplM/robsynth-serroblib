% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:46:35
% EndTime: 2019-02-26 21:46:35
% DurationCPUTime: 0.12s
% Computational Cost: add. (554->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
t71 = qJ(2) + pkin(10) + qJ(4);
t70 = cos(t71);
t69 = sin(t71);
t73 = sin(qJ(1));
t81 = t73 * t69;
t62 = atan2(-t81, -t70);
t58 = sin(t62);
t59 = cos(t62);
t55 = -t58 * t81 - t59 * t70;
t54 = 0.1e1 / t55 ^ 2;
t75 = cos(qJ(1));
t87 = t54 * t75 ^ 2;
t86 = t58 * t70;
t74 = cos(qJ(5));
t77 = t75 * t74;
t72 = sin(qJ(5));
t80 = t73 * t72;
t64 = t70 * t77 + t80;
t61 = 0.1e1 / t64 ^ 2;
t78 = t75 * t72;
t79 = t73 * t74;
t63 = t70 * t78 - t79;
t85 = t61 * t63;
t66 = t69 ^ 2;
t84 = t66 / t70 ^ 2;
t83 = t69 * t75;
t65 = 0.1e1 / (t73 ^ 2 * t84 + 0.1e1);
t82 = t73 * t65;
t76 = t63 ^ 2 * t61 + 0.1e1;
t67 = 0.1e1 / t70;
t60 = 0.1e1 / t64;
t57 = 0.1e1 / t76;
t56 = (0.1e1 + t84) * t82;
t53 = 0.1e1 / t55;
t52 = 0.1e1 / (t66 * t87 + 0.1e1);
t51 = (-t60 * t72 + t74 * t85) * t57 * t83;
t50 = (t70 * t53 - (-t73 * t86 + t59 * t69 + (-t59 * t81 + t86) * t56) * t69 * t54) * t75 * t52;
t1 = [t67 * t65 * t83, t56, 0, t56, 0, 0; (-t53 * t81 - (-t59 * t66 * t67 * t82 + (t65 - 0.1e1) * t69 * t58) * t69 * t87) * t52, t50, 0, t50, 0, 0; ((-t70 * t80 - t77) * t60 - (-t70 * t79 + t78) * t85) * t57, t51, 0, t51, t76 * t57, 0;];
Ja_rot  = t1;
