% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:49:13
% EndTime: 2019-02-26 20:49:13
% DurationCPUTime: 0.11s
% Computational Cost: add. (713->22), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->39)
t74 = qJ(3) + pkin(11) + qJ(5);
t71 = cos(t74);
t70 = sin(t74);
t75 = qJ(1) + pkin(10);
t72 = sin(t75);
t83 = t72 * t70;
t65 = atan2(-t83, -t71);
t63 = sin(t65);
t64 = cos(t65);
t56 = -t63 * t83 - t64 * t71;
t55 = 0.1e1 / t56 ^ 2;
t73 = cos(t75);
t89 = t55 * t73 ^ 2;
t77 = cos(qJ(6));
t79 = t73 * t77;
t76 = sin(qJ(6));
t82 = t72 * t76;
t62 = t71 * t79 + t82;
t60 = 0.1e1 / t62 ^ 2;
t80 = t73 * t76;
t81 = t72 * t77;
t61 = t71 * t80 - t81;
t88 = t60 * t61;
t87 = t63 * t71;
t67 = t70 ^ 2;
t86 = t67 / t71 ^ 2;
t85 = t70 * t73;
t66 = 0.1e1 / (t72 ^ 2 * t86 + 0.1e1);
t84 = t72 * t66;
t78 = t61 ^ 2 * t60 + 0.1e1;
t68 = 0.1e1 / t71;
t59 = 0.1e1 / t62;
t58 = 0.1e1 / t78;
t57 = (0.1e1 + t86) * t84;
t54 = 0.1e1 / t56;
t53 = 0.1e1 / (t67 * t89 + 0.1e1);
t52 = (-t59 * t76 + t77 * t88) * t58 * t85;
t51 = (t71 * t54 - (-t72 * t87 + t64 * t70 + (-t64 * t83 + t87) * t57) * t70 * t55) * t73 * t53;
t1 = [t68 * t66 * t85, 0, t57, 0, t57, 0; (-t54 * t83 - (-t64 * t67 * t68 * t84 + (t66 - 0.1e1) * t70 * t63) * t70 * t89) * t53, 0, t51, 0, t51, 0; ((-t71 * t82 - t79) * t59 - (-t71 * t81 + t80) * t88) * t58, 0, t52, 0, t52, t78 * t58;];
Ja_rot  = t1;
