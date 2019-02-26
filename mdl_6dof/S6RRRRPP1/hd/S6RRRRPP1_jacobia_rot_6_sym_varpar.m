% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP1
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
% Datum: 2019-02-26 22:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:53
% EndTime: 2019-02-26 22:24:53
% DurationCPUTime: 0.12s
% Computational Cost: add. (860->28), mult. (688->69), div. (136->11), fcn. (1034->9), ass. (0->44)
t85 = qJ(2) + qJ(3);
t81 = sin(t85);
t100 = t81 ^ 2;
t82 = cos(t85);
t83 = qJ(4) + pkin(10);
t80 = cos(t83);
t87 = cos(qJ(1));
t89 = t87 * t80;
t79 = sin(t83);
t86 = sin(qJ(1));
t92 = t86 * t79;
t67 = t82 * t92 + t89;
t94 = t81 * t79;
t63 = atan2(-t67, t94);
t60 = sin(t63);
t61 = cos(t63);
t59 = -t60 * t67 + t61 * t94;
t58 = 0.1e1 / t59 ^ 2;
t90 = t87 * t79;
t91 = t86 * t80;
t70 = t82 * t90 - t91;
t99 = t58 * t70;
t97 = t61 * t67;
t96 = t70 ^ 2 * t58;
t74 = 0.1e1 / t79;
t77 = 0.1e1 / t81;
t95 = t74 * t77;
t93 = t81 * t87;
t71 = t82 * t89 + t92;
t66 = 0.1e1 / t71 ^ 2;
t88 = t87 ^ 2 * t100 * t66;
t78 = 0.1e1 / t100;
t75 = 0.1e1 / t79 ^ 2;
t69 = t82 * t91 - t90;
t65 = 0.1e1 / t71;
t64 = 0.1e1 / (0.1e1 + t88);
t62 = 0.1e1 / (t67 ^ 2 * t78 * t75 + 0.1e1);
t57 = 0.1e1 / t59;
t56 = (t67 * t74 * t78 * t82 + t86) * t62;
t55 = (-t65 * t82 * t87 - t80 * t88) * t64;
t54 = 0.1e1 / (0.1e1 + t96);
t53 = (t67 * t75 * t80 - t69 * t74) * t77 * t62;
t52 = (t56 * t97 * t99 + (-t57 * t93 - (t61 * t82 + (-t56 + t86) * t81 * t60) * t99) * t79) * t54;
t1 = [-t70 * t62 * t95, t56, t56, t53, 0, 0; (-t67 * t57 - (-t60 + (t95 * t97 + t60) * t62) * t96) * t54, t52, t52 (t71 * t57 - (t61 * t81 * t80 - t60 * t69 + (-t60 * t94 - t97) * t53) * t99) * t54, 0, 0; (-t66 * t69 * t87 + t65 * t86) * t81 * t64, t55, t55, -t70 * t66 * t64 * t93, 0, 0;];
Ja_rot  = t1;
