% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR10
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR10_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:17
% EndTime: 2019-02-26 21:59:18
% DurationCPUTime: 0.22s
% Computational Cost: add. (878->38), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
t84 = cos(pkin(6));
t86 = sin(qJ(2));
t90 = cos(qJ(1));
t94 = t90 * t86;
t87 = sin(qJ(1));
t89 = cos(qJ(2));
t95 = t87 * t89;
t75 = t84 * t94 + t95;
t82 = pkin(12) + qJ(4);
t80 = sin(t82);
t81 = cos(t82);
t83 = sin(pkin(6));
t97 = t83 * t90;
t64 = t75 * t80 + t81 * t97;
t100 = t83 * t86;
t72 = t80 * t100 - t84 * t81;
t63 = atan2(-t64, t72);
t56 = sin(t63);
t57 = cos(t63);
t54 = -t56 * t64 + t57 * t72;
t53 = 0.1e1 / t54 ^ 2;
t93 = t90 * t89;
t96 = t87 * t86;
t77 = -t84 * t96 + t93;
t99 = t83 * t87;
t68 = t77 * t80 - t81 * t99;
t107 = t53 * t68;
t106 = t57 * t64;
t76 = t84 * t95 + t94;
t85 = sin(qJ(5));
t102 = t76 * t85;
t69 = t77 * t81 + t80 * t99;
t88 = cos(qJ(5));
t62 = t69 * t88 + t102;
t59 = 0.1e1 / t62 ^ 2;
t101 = t76 * t88;
t61 = t69 * t85 - t101;
t105 = t59 * t61;
t71 = 0.1e1 / t72 ^ 2;
t104 = t64 * t71;
t103 = t68 ^ 2 * t53;
t98 = t83 * t89;
t92 = t61 ^ 2 * t59 + 0.1e1;
t66 = t75 * t81 - t80 * t97;
t91 = -t56 * t72 - t106;
t74 = t84 * t93 - t96;
t73 = t81 * t100 + t84 * t80;
t70 = 0.1e1 / t72;
t60 = 0.1e1 / (t64 ^ 2 * t71 + 0.1e1);
t58 = 0.1e1 / t62;
t55 = 0.1e1 / t92;
t52 = 0.1e1 / t54;
t51 = 0.1e1 / (0.1e1 + t103);
t50 = (t98 * t104 - t70 * t74) * t80 * t60;
t49 = (t73 * t104 - t66 * t70) * t60;
t1 = [-t68 * t70 * t60, t50, 0, t49, 0, 0; (-t64 * t52 - (-t56 + (t70 * t106 + t56) * t60) * t103) * t51 (-t76 * t80 * t52 - ((-t56 * t74 + t57 * t98) * t80 + t91 * t50) * t107) * t51, 0 (t69 * t52 - (t91 * t49 - t56 * t66 + t57 * t73) * t107) * t51, 0, 0; ((-t66 * t85 - t74 * t88) * t58 - (-t66 * t88 + t74 * t85) * t105) * t55 ((-t81 * t102 - t77 * t88) * t58 - (-t81 * t101 + t77 * t85) * t105) * t55, 0 (t88 * t105 - t85 * t58) * t68 * t55, t92 * t55, 0;];
Ja_rot  = t1;
