% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_rot = S6RRPRRP2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:46:35
% EndTime: 2019-02-26 21:46:35
% DurationCPUTime: 0.14s
% Computational Cost: add. (740->27), mult. (688->69), div. (136->11), fcn. (1034->9), ass. (0->43)
t82 = qJ(2) + pkin(10) + qJ(4);
t78 = sin(t82);
t102 = t78 ^ 2;
t79 = cos(t82);
t88 = cos(qJ(5));
t89 = cos(qJ(1));
t91 = t89 * t88;
t86 = sin(qJ(5));
t87 = sin(qJ(1));
t94 = t87 * t86;
t70 = t79 * t94 + t91;
t96 = t78 * t86;
t67 = atan2(-t70, t96);
t63 = sin(t67);
t64 = cos(t67);
t62 = -t63 * t70 + t64 * t96;
t61 = 0.1e1 / t62 ^ 2;
t92 = t89 * t86;
t93 = t87 * t88;
t73 = t79 * t92 - t93;
t101 = t61 * t73;
t99 = t64 * t70;
t98 = t73 ^ 2 * t61;
t76 = 0.1e1 / t78;
t83 = 0.1e1 / t86;
t97 = t76 * t83;
t95 = t78 * t89;
t74 = t79 * t91 + t94;
t69 = 0.1e1 / t74 ^ 2;
t90 = t89 ^ 2 * t102 * t69;
t84 = 0.1e1 / t86 ^ 2;
t77 = 0.1e1 / t102;
t72 = t79 * t93 - t92;
t68 = 0.1e1 / t74;
t66 = 0.1e1 / (t70 ^ 2 * t77 * t84 + 0.1e1);
t65 = 0.1e1 / (0.1e1 + t90);
t60 = 0.1e1 / t62;
t59 = (t70 * t77 * t79 * t83 + t87) * t66;
t58 = 0.1e1 / (0.1e1 + t98);
t57 = (t70 * t84 * t88 - t72 * t83) * t76 * t66;
t56 = (-t68 * t79 * t89 - t88 * t90) * t65;
t55 = (t59 * t99 * t101 + (-t60 * t95 - (t64 * t79 + (-t59 + t87) * t78 * t63) * t101) * t86) * t58;
t1 = [-t73 * t66 * t97, t59, 0, t59, t57, 0; (-t70 * t60 - (-t63 + (t97 * t99 + t63) * t66) * t98) * t58, t55, 0, t55 (t74 * t60 - (t64 * t78 * t88 - t63 * t72 + (-t63 * t96 - t99) * t57) * t101) * t58, 0; (-t69 * t72 * t89 + t68 * t87) * t78 * t65, t56, 0, t56, -t73 * t69 * t65 * t95, 0;];
Ja_rot  = t1;
