% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:26:59
% EndTime: 2019-02-26 22:26:59
% DurationCPUTime: 0.13s
% Computational Cost: add. (1004->28), mult. (689->68), div. (139->11), fcn. (1046->9), ass. (0->42)
t81 = sin(qJ(2));
t96 = t81 ^ 2;
t76 = qJ(3) + qJ(4) + pkin(10);
t74 = sin(t76);
t75 = cos(t76);
t84 = cos(qJ(1));
t86 = t84 * t75;
t82 = sin(qJ(1));
t83 = cos(qJ(2));
t88 = t82 * t83;
t65 = t74 * t88 + t86;
t90 = t81 * t74;
t61 = atan2(-t65, t90);
t58 = sin(t61);
t59 = cos(t61);
t56 = -t58 * t65 + t59 * t90;
t55 = 0.1e1 / t56 ^ 2;
t87 = t84 * t74;
t68 = -t82 * t75 + t83 * t87;
t95 = t55 * t68;
t93 = t59 * t65;
t92 = t68 ^ 2 * t55;
t72 = 0.1e1 / t74;
t78 = 0.1e1 / t81;
t91 = t72 * t78;
t89 = t81 * t84;
t69 = t82 * t74 + t83 * t86;
t64 = 0.1e1 / t69 ^ 2;
t85 = t84 ^ 2 * t96 * t64;
t79 = 0.1e1 / t96;
t73 = 0.1e1 / t74 ^ 2;
t67 = t75 * t88 - t87;
t63 = 0.1e1 / t69;
t62 = 0.1e1 / (0.1e1 + t85);
t60 = 0.1e1 / (t65 ^ 2 * t79 * t73 + 0.1e1);
t57 = t68 * t64 * t62 * t89;
t54 = 0.1e1 / t56;
t53 = (t65 * t72 * t79 * t83 + t82) * t60;
t52 = 0.1e1 / (0.1e1 + t92);
t51 = (t65 * t73 * t75 - t67 * t72) * t78 * t60;
t50 = (t69 * t54 - (t59 * t81 * t75 - t58 * t67 + (-t58 * t90 - t93) * t51) * t95) * t52;
t1 = [-t68 * t60 * t91, t53, t51, t51, 0, 0; (-t65 * t54 - (-t58 + (t91 * t93 + t58) * t60) * t92) * t52 (t53 * t93 * t95 + (-t54 * t89 - (t59 * t83 + (-t53 + t82) * t81 * t58) * t95) * t74) * t52, t50, t50, 0, 0; (-t64 * t67 * t84 + t63 * t82) * t81 * t62 (-t63 * t83 * t84 - t75 * t85) * t62, -t57, -t57, 0, 0;];
Ja_rot  = t1;
