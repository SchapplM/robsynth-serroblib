% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:11:19
% EndTime: 2019-02-26 21:11:19
% DurationCPUTime: 0.15s
% Computational Cost: add. (873->29), mult. (689->67), div. (139->11), fcn. (1046->9), ass. (0->44)
t80 = pkin(10) + qJ(3);
t76 = sin(t80);
t97 = t76 ^ 2;
t77 = cos(t80);
t82 = qJ(4) + qJ(5);
t78 = sin(t82);
t83 = sin(qJ(1));
t88 = t83 * t78;
t79 = cos(t82);
t84 = cos(qJ(1));
t89 = t79 * t84;
t64 = t77 * t88 + t89;
t91 = t76 * t78;
t60 = atan2(-t64, t91);
t57 = sin(t60);
t58 = cos(t60);
t55 = -t57 * t64 + t58 * t91;
t54 = 0.1e1 / t55 ^ 2;
t86 = t84 * t78;
t87 = t83 * t79;
t67 = t77 * t86 - t87;
t96 = t54 * t67;
t95 = t54 * t67 ^ 2;
t93 = t58 * t64;
t72 = 0.1e1 / t76;
t74 = 0.1e1 / t78;
t92 = t72 * t74;
t90 = t76 * t84;
t68 = t77 * t89 + t88;
t63 = 0.1e1 / t68 ^ 2;
t85 = t63 * t97 * t84 ^ 2;
t75 = 0.1e1 / t78 ^ 2;
t73 = 0.1e1 / t97;
t66 = t77 * t87 - t86;
t62 = 0.1e1 / t68;
t61 = 0.1e1 / (0.1e1 + t85);
t59 = 0.1e1 / (t64 ^ 2 * t73 * t75 + 0.1e1);
t56 = t67 * t63 * t61 * t90;
t53 = 0.1e1 / t55;
t52 = (t64 * t73 * t74 * t77 + t83) * t59;
t51 = 0.1e1 / (0.1e1 + t95);
t50 = (t64 * t75 * t79 - t66 * t74) * t72 * t59;
t49 = (t68 * t53 - (t58 * t76 * t79 - t57 * t66 + (-t57 * t91 - t93) * t50) * t96) * t51;
t1 = [-t67 * t59 * t92, 0, t52, t50, t50, 0; (-t64 * t53 - (-t57 + (t92 * t93 + t57) * t59) * t95) * t51, 0 (t52 * t93 * t96 + (-t53 * t90 - (t58 * t77 + (-t52 + t83) * t76 * t57) * t96) * t78) * t51, t49, t49, 0; (-t63 * t66 * t84 + t62 * t83) * t76 * t61, 0 (-t62 * t77 * t84 - t79 * t85) * t61, -t56, -t56, 0;];
Ja_rot  = t1;
