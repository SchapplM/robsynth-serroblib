% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:10:35
% EndTime: 2019-02-26 22:10:36
% DurationCPUTime: 0.15s
% Computational Cost: add. (860->28), mult. (688->69), div. (136->11), fcn. (1034->9), ass. (0->44)
t82 = qJ(2) + qJ(3);
t78 = sin(t82);
t97 = t78 ^ 2;
t80 = pkin(10) + qJ(5);
t77 = cos(t80);
t79 = cos(t82);
t84 = cos(qJ(1));
t76 = sin(t80);
t83 = sin(qJ(1));
t88 = t83 * t76;
t64 = t77 * t84 + t79 * t88;
t91 = t78 * t76;
t60 = atan2(-t64, t91);
t57 = sin(t60);
t58 = cos(t60);
t56 = -t57 * t64 + t58 * t91;
t55 = 0.1e1 / t56 ^ 2;
t86 = t84 * t76;
t87 = t83 * t77;
t67 = t79 * t86 - t87;
t96 = t55 * t67;
t95 = t55 * t67 ^ 2;
t93 = t58 * t64;
t71 = 0.1e1 / t76;
t74 = 0.1e1 / t78;
t92 = t71 * t74;
t90 = t78 * t84;
t89 = t79 * t84;
t68 = t77 * t89 + t88;
t63 = 0.1e1 / t68 ^ 2;
t85 = t63 * t97 * t84 ^ 2;
t75 = 0.1e1 / t97;
t72 = 0.1e1 / t76 ^ 2;
t66 = t79 * t87 - t86;
t62 = 0.1e1 / t68;
t61 = 0.1e1 / (0.1e1 + t85);
t59 = 0.1e1 / (t64 ^ 2 * t72 * t75 + 0.1e1);
t54 = 0.1e1 / t56;
t53 = (t64 * t71 * t75 * t79 + t83) * t59;
t52 = (-t62 * t89 - t77 * t85) * t61;
t51 = 0.1e1 / (0.1e1 + t95);
t50 = (t64 * t72 * t77 - t66 * t71) * t74 * t59;
t49 = (t53 * t93 * t96 + (-t54 * t90 - (t58 * t79 + (-t53 + t83) * t78 * t57) * t96) * t76) * t51;
t1 = [-t67 * t59 * t92, t53, t53, 0, t50, 0; (-t64 * t54 - (-t57 + (t92 * t93 + t57) * t59) * t95) * t51, t49, t49, 0 (t68 * t54 - (t58 * t77 * t78 - t57 * t66 + (-t57 * t91 - t93) * t50) * t96) * t51, 0; (-t63 * t66 * t84 + t62 * t83) * t78 * t61, t52, t52, 0, -t67 * t63 * t61 * t90, 0;];
Ja_rot  = t1;
