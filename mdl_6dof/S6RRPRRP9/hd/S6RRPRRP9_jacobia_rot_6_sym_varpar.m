% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:50:28
% EndTime: 2019-02-26 21:50:28
% DurationCPUTime: 0.19s
% Computational Cost: add. (878->38), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
t86 = cos(pkin(6));
t88 = sin(qJ(2));
t92 = cos(qJ(1));
t96 = t92 * t88;
t89 = sin(qJ(1));
t91 = cos(qJ(2));
t97 = t89 * t91;
t77 = t86 * t96 + t97;
t84 = pkin(11) + qJ(4);
t82 = sin(t84);
t83 = cos(t84);
t85 = sin(pkin(6));
t99 = t85 * t92;
t66 = t77 * t82 + t83 * t99;
t102 = t85 * t88;
t74 = t102 * t82 - t83 * t86;
t65 = atan2(-t66, t74);
t58 = sin(t65);
t59 = cos(t65);
t56 = -t58 * t66 + t59 * t74;
t55 = 0.1e1 / t56 ^ 2;
t101 = t85 * t89;
t95 = t92 * t91;
t98 = t89 * t88;
t79 = -t86 * t98 + t95;
t70 = -t101 * t83 + t79 * t82;
t109 = t55 * t70;
t108 = t59 * t66;
t78 = t86 * t97 + t96;
t87 = sin(qJ(5));
t104 = t78 * t87;
t71 = t101 * t82 + t79 * t83;
t90 = cos(qJ(5));
t64 = t71 * t90 + t104;
t61 = 0.1e1 / t64 ^ 2;
t103 = t78 * t90;
t63 = t71 * t87 - t103;
t107 = t61 * t63;
t73 = 0.1e1 / t74 ^ 2;
t106 = t66 * t73;
t105 = t70 ^ 2 * t55;
t100 = t85 * t91;
t94 = t61 * t63 ^ 2 + 0.1e1;
t68 = t77 * t83 - t82 * t99;
t93 = -t58 * t74 - t108;
t76 = t86 * t95 - t98;
t75 = t102 * t83 + t82 * t86;
t72 = 0.1e1 / t74;
t62 = 0.1e1 / (t66 ^ 2 * t73 + 0.1e1);
t60 = 0.1e1 / t64;
t57 = 0.1e1 / t94;
t54 = 0.1e1 / t56;
t53 = 0.1e1 / (0.1e1 + t105);
t52 = (t100 * t106 - t72 * t76) * t82 * t62;
t51 = (t106 * t75 - t68 * t72) * t62;
t1 = [-t70 * t72 * t62, t52, 0, t51, 0, 0; (-t66 * t54 - (-t58 + (t108 * t72 + t58) * t62) * t105) * t53 (-t78 * t82 * t54 - ((t100 * t59 - t58 * t76) * t82 + t93 * t52) * t109) * t53, 0 (t71 * t54 - (t51 * t93 - t58 * t68 + t59 * t75) * t109) * t53, 0, 0; ((-t68 * t87 - t76 * t90) * t60 - (-t68 * t90 + t76 * t87) * t107) * t57 ((-t104 * t83 - t79 * t90) * t60 - (-t103 * t83 + t79 * t87) * t107) * t57, 0 (t107 * t90 - t60 * t87) * t70 * t57, t94 * t57, 0;];
Ja_rot  = t1;
