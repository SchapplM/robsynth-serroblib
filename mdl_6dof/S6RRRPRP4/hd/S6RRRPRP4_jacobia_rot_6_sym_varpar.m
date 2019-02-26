% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:11:05
% EndTime: 2019-02-26 22:11:05
% DurationCPUTime: 0.13s
% Computational Cost: add. (477->26), mult. (688->69), div. (136->11), fcn. (1034->9), ass. (0->43)
t80 = qJ(2) + qJ(3);
t76 = cos(t80);
t97 = t76 ^ 2;
t75 = sin(t80);
t81 = sin(qJ(5));
t84 = cos(qJ(1));
t87 = t84 * t81;
t82 = sin(qJ(1));
t83 = cos(qJ(5));
t88 = t82 * t83;
t69 = t75 * t88 + t87;
t91 = t76 * t83;
t64 = atan2(t69, t91);
t60 = sin(t64);
t61 = cos(t64);
t59 = t60 * t69 + t61 * t91;
t58 = 0.1e1 / t59 ^ 2;
t86 = t84 * t83;
t89 = t82 * t81;
t67 = -t75 * t86 + t89;
t96 = t58 * t67;
t94 = t61 * t69;
t93 = t67 ^ 2 * t58;
t73 = 0.1e1 / t76;
t77 = 0.1e1 / t83;
t92 = t73 * t77;
t90 = t76 * t84;
t68 = t75 * t87 + t88;
t66 = 0.1e1 / t68 ^ 2;
t85 = t84 ^ 2 * t97 * t66;
t78 = 0.1e1 / t83 ^ 2;
t74 = 0.1e1 / t97;
t70 = -t75 * t89 + t86;
t65 = 0.1e1 / t68;
t63 = 0.1e1 / (t69 ^ 2 * t74 * t78 + 0.1e1);
t62 = 0.1e1 / (0.1e1 + t85);
t57 = 0.1e1 / t59;
t56 = (t69 * t74 * t75 * t77 + t82) * t63;
t55 = 0.1e1 / (0.1e1 + t93);
t54 = (t69 * t78 * t81 + t70 * t77) * t73 * t63;
t53 = (t65 * t75 * t84 + t81 * t85) * t62;
t52 = (-t56 * t94 * t96 + (-t57 * t90 - (-t61 * t75 + (-t56 + t82) * t60 * t76) * t96) * t83) * t55;
t1 = [-t67 * t63 * t92, t56, t56, 0, t54, 0; (t69 * t57 - (-t60 + (-t92 * t94 + t60) * t63) * t93) * t55, t52, t52, 0 (t68 * t57 - (-t61 * t76 * t81 + t60 * t70 + (-t60 * t91 + t94) * t54) * t96) * t55, 0; (t66 * t70 * t84 + t65 * t82) * t76 * t62, t53, t53, 0, -t67 * t66 * t62 * t90, 0;];
Ja_rot  = t1;
