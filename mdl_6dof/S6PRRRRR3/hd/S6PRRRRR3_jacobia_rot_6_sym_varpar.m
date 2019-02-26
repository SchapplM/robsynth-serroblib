% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_jacobia_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:19:52
% EndTime: 2019-02-26 20:19:52
% DurationCPUTime: 0.12s
% Computational Cost: add. (548->31), mult. (1115->79), div. (75->9), fcn. (1573->13), ass. (0->51)
t85 = sin(pkin(12));
t87 = cos(pkin(12));
t92 = cos(qJ(2));
t88 = cos(pkin(6));
t90 = sin(qJ(2));
t96 = t88 * t90;
t76 = t85 * t92 + t87 * t96;
t89 = sin(qJ(3));
t86 = sin(pkin(6));
t91 = cos(qJ(3));
t98 = t86 * t91;
t68 = t76 * t89 + t87 * t98;
t99 = t86 * t89;
t79 = -t88 * t91 + t90 * t99;
t67 = atan2(-t68, t79);
t64 = sin(t67);
t65 = cos(t67);
t58 = -t64 * t68 + t65 * t79;
t57 = 0.1e1 / t58 ^ 2;
t78 = -t85 * t96 + t87 * t92;
t71 = t78 * t89 - t85 * t98;
t103 = t57 * t71;
t72 = t78 * t91 + t85 * t99;
t95 = t88 * t92;
t77 = t85 * t95 + t87 * t90;
t84 = qJ(4) + qJ(5) + qJ(6);
t82 = sin(t84);
t83 = cos(t84);
t63 = t72 * t83 + t77 * t82;
t61 = 0.1e1 / t63 ^ 2;
t62 = t72 * t82 - t77 * t83;
t102 = t61 * t62;
t74 = 0.1e1 / t79 ^ 2;
t101 = t68 * t74;
t100 = t77 * t91;
t97 = t86 * t92;
t94 = t62 ^ 2 * t61 + 0.1e1;
t93 = -t64 * t79 - t65 * t68;
t80 = t88 * t89 + t90 * t98;
t75 = -t85 * t90 + t87 * t95;
t73 = 0.1e1 / t79;
t70 = t76 * t91 - t87 * t99;
t66 = 0.1e1 / (t68 ^ 2 * t74 + 0.1e1);
t60 = 0.1e1 / t63;
t59 = 0.1e1 / t94;
t56 = 0.1e1 / t58;
t55 = 0.1e1 / (t71 ^ 2 * t57 + 0.1e1);
t54 = (t97 * t101 - t73 * t75) * t89 * t66;
t53 = (t80 * t101 - t70 * t73) * t66;
t52 = t94 * t59;
t1 = [0, t54, t53, 0, 0, 0; 0 (-t77 * t89 * t56 - ((-t64 * t75 + t65 * t97) * t89 + t93 * t54) * t103) * t55 (t72 * t56 - (t93 * t53 - t64 * t70 + t65 * t80) * t103) * t55, 0, 0, 0; 0 ((-t82 * t100 - t78 * t83) * t60 - (-t83 * t100 + t78 * t82) * t102) * t59 (t83 * t102 - t60 * t82) * t71 * t59, t52, t52, t52;];
Ja_rot  = t1;
