% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPP3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:01
% EndTime: 2019-02-26 20:10:01
% DurationCPUTime: 0.18s
% Computational Cost: add. (707->41), mult. (1976->104), div. (87->9), fcn. (2786->13), ass. (0->56)
t89 = sin(pkin(6));
t93 = sin(qJ(3));
t104 = t89 * t93;
t91 = cos(pkin(6));
t94 = sin(qJ(2));
t102 = t91 * t94;
t88 = sin(pkin(10));
t90 = cos(pkin(10));
t97 = cos(qJ(2));
t82 = t90 * t102 + t88 * t97;
t96 = cos(qJ(3));
t73 = -t90 * t104 + t82 * t96;
t101 = t91 * t97;
t81 = -t90 * t101 + t88 * t94;
t92 = sin(qJ(4));
t95 = cos(qJ(4));
t65 = t73 * t95 + t81 * t92;
t103 = t89 * t96;
t86 = t94 * t103 + t91 * t93;
t79 = -t89 * t97 * t92 + t86 * t95;
t63 = atan2(-t65, t79);
t60 = sin(t63);
t61 = cos(t63);
t58 = -t60 * t65 + t61 * t79;
t57 = 0.1e1 / t58 ^ 2;
t83 = t88 * t101 + t90 * t94;
t105 = t83 * t92;
t84 = -t88 * t102 + t90 * t97;
t75 = t88 * t104 + t84 * t96;
t68 = t75 * t95 + t105;
t108 = t57 * t68;
t77 = 0.1e1 / t79 ^ 2;
t107 = t65 * t77;
t67 = -t75 * t92 + t83 * t95;
t74 = -t88 * t103 + t84 * t93;
t71 = 0.1e1 / t74 ^ 2;
t106 = t67 * t71;
t100 = t95 * t96;
t99 = t95 * t97;
t98 = -t60 * t79 - t61 * t65;
t85 = -t94 * t104 + t91 * t96;
t80 = (t92 * t94 + t96 * t99) * t89;
t78 = -t86 * t92 - t89 * t99;
t76 = 0.1e1 / t79;
t72 = -t90 * t103 - t82 * t93;
t70 = 0.1e1 / t74;
t69 = -t81 * t100 + t82 * t92;
t64 = t73 * t92 - t81 * t95;
t62 = 0.1e1 / (t65 ^ 2 * t77 + 0.1e1);
t59 = 0.1e1 / (t67 ^ 2 * t71 + 0.1e1);
t56 = 0.1e1 / t58;
t55 = 0.1e1 / (t57 * t68 ^ 2 + 0.1e1);
t54 = (t85 * t107 - t72 * t76) * t95 * t62;
t53 = (t80 * t107 - t69 * t76) * t62;
t52 = (t78 * t107 + t64 * t76) * t62;
t1 = [0, t53, t54, t52, 0, 0; 0 ((-t83 * t100 + t84 * t92) * t56 - (t98 * t53 - t60 * t69 + t61 * t80) * t108) * t55 (-t74 * t95 * t56 - ((-t60 * t72 + t61 * t85) * t95 + t98 * t54) * t108) * t55 (t67 * t56 - (t98 * t52 + t60 * t64 + t61 * t78) * t108) * t55, 0, 0; 0 ((t96 * t105 + t84 * t95) * t70 + t83 * t93 * t106) * t59 (t70 * t74 * t92 - t75 * t106) * t59, -t68 * t70 * t59, 0, 0;];
Ja_rot  = t1;
