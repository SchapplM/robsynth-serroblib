% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:13
% EndTime: 2019-02-26 21:03:13
% DurationCPUTime: 0.08s
% Computational Cost: add. (520->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
t71 = pkin(10) + qJ(3) + qJ(4);
t69 = sin(t71);
t70 = cos(t71);
t73 = sin(qJ(1));
t81 = t73 * t70;
t64 = atan2(-t81, t69);
t58 = sin(t64);
t59 = cos(t64);
t55 = -t58 * t81 + t59 * t69;
t54 = 0.1e1 / t55 ^ 2;
t75 = cos(qJ(1));
t87 = t54 * t75 ^ 2;
t86 = t58 * t69;
t72 = sin(qJ(6));
t78 = t75 * t72;
t74 = cos(qJ(6));
t79 = t73 * t74;
t63 = t69 * t78 + t79;
t61 = 0.1e1 / t63 ^ 2;
t77 = t75 * t74;
t80 = t73 * t72;
t62 = -t69 * t77 + t80;
t85 = t61 * t62;
t68 = t70 ^ 2;
t84 = 0.1e1 / t69 ^ 2 * t68;
t83 = t70 * t75;
t65 = 0.1e1 / (t73 ^ 2 * t84 + 0.1e1);
t82 = t73 * t65;
t76 = t62 ^ 2 * t61 + 0.1e1;
t66 = 0.1e1 / t69;
t60 = 0.1e1 / t63;
t57 = 0.1e1 / t76;
t56 = (0.1e1 + t84) * t82;
t53 = 0.1e1 / t55;
t52 = 0.1e1 / (t68 * t87 + 0.1e1);
t51 = (-t60 * t74 - t72 * t85) * t57 * t83;
t50 = (-t69 * t53 - (t73 * t86 + t59 * t70 + (-t59 * t81 - t86) * t56) * t70 * t54) * t75 * t52;
t1 = [-t66 * t65 * t83, 0, t56, t56, 0, 0; (-t53 * t81 - (t59 * t66 * t68 * t82 + (t65 - 0.1e1) * t70 * t58) * t70 * t87) * t52, 0, t50, t50, 0, 0; ((t69 * t79 + t78) * t60 - (-t69 * t80 + t77) * t85) * t57, 0, t51, t51, 0, t76 * t57;];
Ja_rot  = t1;
