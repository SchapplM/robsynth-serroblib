% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPPR4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:00
% EndTime: 2019-02-26 20:00:00
% DurationCPUTime: 0.15s
% Computational Cost: add. (475->34), mult. (1334->86), div. (60->9), fcn. (1876->13), ass. (0->50)
t77 = sin(pkin(10));
t80 = cos(pkin(10));
t85 = cos(qJ(2));
t81 = cos(pkin(6));
t83 = sin(qJ(2));
t89 = t81 * t83;
t72 = t77 * t85 + t80 * t89;
t76 = sin(pkin(11));
t79 = cos(pkin(11));
t84 = cos(qJ(3));
t88 = t81 * t85;
t86 = -t77 * t83 + t80 * t88;
t78 = sin(pkin(6));
t82 = sin(qJ(3));
t91 = t78 * t82;
t58 = (t72 * t84 - t80 * t91) * t76 + t86 * t79;
t90 = t78 * t84;
t66 = (t81 * t82 + t83 * t90) * t76 + t78 * t85 * t79;
t55 = atan2(-t58, t66);
t51 = sin(t55);
t52 = cos(t55);
t50 = -t51 * t58 + t52 * t66;
t49 = 0.1e1 / t50 ^ 2;
t74 = -t77 * t89 + t80 * t85;
t69 = t74 * t84 + t77 * t91;
t73 = t77 * t88 + t80 * t83;
t93 = t73 * t79;
t60 = t69 * t76 - t93;
t96 = t49 * t60;
t64 = 0.1e1 / t66 ^ 2;
t95 = t58 * t64;
t61 = t69 * t79 + t73 * t76;
t57 = 0.1e1 / t61 ^ 2;
t68 = -t74 * t82 + t77 * t90;
t94 = t68 ^ 2 * t57;
t92 = t76 * t84;
t87 = -t51 * t66 - t52 * t58;
t75 = t81 * t84 - t83 * t91;
t70 = (-t79 * t83 + t85 * t92) * t78;
t67 = -t72 * t82 - t80 * t90;
t63 = 0.1e1 / t66;
t62 = -t72 * t79 + t86 * t92;
t56 = 0.1e1 / t61;
t54 = 0.1e1 / (t58 ^ 2 * t64 + 0.1e1);
t53 = 0.1e1 / (0.1e1 + t94);
t48 = 0.1e1 / t50;
t47 = 0.1e1 / (t60 ^ 2 * t49 + 0.1e1);
t46 = (-t63 * t67 + t75 * t95) * t76 * t54;
t45 = (-t62 * t63 + t70 * t95) * t54;
t1 = [0, t45, t46, 0, 0, 0; 0 ((-t73 * t92 - t74 * t79) * t48 - (t45 * t87 - t51 * t62 + t52 * t70) * t96) * t47 (t68 * t76 * t48 - ((-t51 * t67 + t52 * t75) * t76 + t87 * t46) * t96) * t47, 0, 0, 0; 0 (t73 * t82 * t56 - (t74 * t76 - t84 * t93) * t68 * t57) * t53 (-t56 * t69 - t79 * t94) * t53, 0, 0, 0;];
Ja_rot  = t1;
