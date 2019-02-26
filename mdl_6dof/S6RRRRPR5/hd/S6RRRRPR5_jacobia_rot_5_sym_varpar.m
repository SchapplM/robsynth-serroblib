% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:32:53
% EndTime: 2019-02-26 22:32:53
% DurationCPUTime: 0.15s
% Computational Cost: add. (477->27), mult. (688->69), div. (136->11), fcn. (1034->9), ass. (0->43)
t82 = qJ(2) + qJ(3);
t77 = sin(t82);
t99 = t77 ^ 2;
t78 = cos(t82);
t85 = cos(qJ(4));
t86 = cos(qJ(1));
t88 = t86 * t85;
t83 = sin(qJ(4));
t84 = sin(qJ(1));
t91 = t84 * t83;
t67 = t78 * t91 + t88;
t93 = t77 * t83;
t64 = atan2(-t67, t93);
t60 = sin(t64);
t61 = cos(t64);
t59 = -t60 * t67 + t61 * t93;
t58 = 0.1e1 / t59 ^ 2;
t89 = t86 * t83;
t90 = t84 * t85;
t70 = t78 * t89 - t90;
t98 = t58 * t70;
t96 = t61 * t67;
t95 = t70 ^ 2 * t58;
t75 = 0.1e1 / t77;
t79 = 0.1e1 / t83;
t94 = t75 * t79;
t92 = t77 * t86;
t71 = t78 * t88 + t91;
t66 = 0.1e1 / t71 ^ 2;
t87 = t86 ^ 2 * t99 * t66;
t80 = 0.1e1 / t83 ^ 2;
t76 = 0.1e1 / t99;
t69 = t78 * t90 - t89;
t65 = 0.1e1 / t71;
t63 = 0.1e1 / (t67 ^ 2 * t76 * t80 + 0.1e1);
t62 = 0.1e1 / (0.1e1 + t87);
t57 = 0.1e1 / t59;
t56 = (t67 * t76 * t78 * t79 + t84) * t63;
t55 = 0.1e1 / (0.1e1 + t95);
t54 = (t67 * t80 * t85 - t69 * t79) * t75 * t63;
t53 = (-t65 * t78 * t86 - t85 * t87) * t62;
t52 = (t56 * t96 * t98 + (-t57 * t92 - (t61 * t78 + (-t56 + t84) * t60 * t77) * t98) * t83) * t55;
t1 = [-t70 * t63 * t94, t56, t56, t54, 0, 0; (-t67 * t57 - (-t60 + (t94 * t96 + t60) * t63) * t95) * t55, t52, t52 (t71 * t57 - (t61 * t77 * t85 - t60 * t69 + (-t60 * t93 - t96) * t54) * t98) * t55, 0, 0; (-t66 * t69 * t86 + t65 * t84) * t77 * t62, t53, t53, -t70 * t66 * t62 * t92, 0, 0;];
Ja_rot  = t1;
