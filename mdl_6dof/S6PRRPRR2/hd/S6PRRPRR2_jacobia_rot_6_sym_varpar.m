% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:04:42
% EndTime: 2019-02-26 20:04:42
% DurationCPUTime: 0.20s
% Computational Cost: add. (748->32), mult. (1032->79), div. (70->9), fcn. (1461->13), ass. (0->54)
t90 = sin(pkin(6));
t91 = cos(pkin(11));
t101 = t90 * t91;
t89 = sin(pkin(11));
t94 = cos(qJ(2));
t92 = cos(pkin(6));
t93 = sin(qJ(2));
t98 = t92 * t93;
t79 = t89 * t94 + t91 * t98;
t87 = qJ(3) + pkin(12);
t83 = sin(t87);
t84 = cos(t87);
t69 = t84 * t101 + t79 * t83;
t100 = t90 * t93;
t76 = t83 * t100 - t92 * t84;
t68 = atan2(-t69, t76);
t65 = sin(t68);
t66 = cos(t68);
t59 = -t65 * t69 + t66 * t76;
t58 = 0.1e1 / t59 ^ 2;
t102 = t89 * t90;
t81 = -t89 * t98 + t91 * t94;
t72 = -t84 * t102 + t81 * t83;
t107 = t58 * t72;
t97 = t92 * t94;
t80 = t89 * t97 + t91 * t93;
t88 = qJ(5) + qJ(6);
t85 = sin(t88);
t104 = t80 * t85;
t73 = t83 * t102 + t81 * t84;
t86 = cos(t88);
t64 = t73 * t86 + t104;
t62 = 0.1e1 / t64 ^ 2;
t103 = t80 * t86;
t63 = t73 * t85 - t103;
t106 = t62 * t63;
t75 = 0.1e1 / t76 ^ 2;
t105 = t69 * t75;
t99 = t90 * t94;
t96 = t63 ^ 2 * t62 + 0.1e1;
t95 = -t65 * t76 - t66 * t69;
t78 = -t89 * t93 + t91 * t97;
t77 = t84 * t100 + t92 * t83;
t74 = 0.1e1 / t76;
t71 = -t83 * t101 + t79 * t84;
t67 = 0.1e1 / (t69 ^ 2 * t75 + 0.1e1);
t61 = 0.1e1 / t64;
t60 = 0.1e1 / t96;
t57 = 0.1e1 / t59;
t56 = 0.1e1 / (t72 ^ 2 * t58 + 0.1e1);
t55 = (t99 * t105 - t74 * t78) * t83 * t67;
t54 = (t77 * t105 - t71 * t74) * t67;
t53 = t96 * t60;
t1 = [0, t55, t54, 0, 0, 0; 0 (-t80 * t83 * t57 - ((-t65 * t78 + t66 * t99) * t83 + t95 * t55) * t107) * t56 (t73 * t57 - (t95 * t54 - t65 * t71 + t66 * t77) * t107) * t56, 0, 0, 0; 0 ((-t84 * t104 - t81 * t86) * t61 - (-t84 * t103 + t81 * t85) * t106) * t60 (t86 * t106 - t61 * t85) * t72 * t60, 0, t53, t53;];
Ja_rot  = t1;
