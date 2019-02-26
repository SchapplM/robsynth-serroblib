% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP4
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
% Datum: 2019-02-26 21:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:47:45
% EndTime: 2019-02-26 21:47:45
% DurationCPUTime: 0.15s
% Computational Cost: add. (873->29), mult. (689->67), div. (139->11), fcn. (1046->9), ass. (0->44)
t83 = qJ(2) + pkin(10);
t79 = sin(t83);
t100 = t79 ^ 2;
t80 = cos(t83);
t85 = qJ(4) + qJ(5);
t81 = sin(t85);
t86 = sin(qJ(1));
t91 = t86 * t81;
t82 = cos(t85);
t87 = cos(qJ(1));
t92 = t82 * t87;
t67 = t80 * t91 + t92;
t94 = t79 * t81;
t63 = atan2(-t67, t94);
t60 = sin(t63);
t61 = cos(t63);
t58 = -t60 * t67 + t61 * t94;
t57 = 0.1e1 / t58 ^ 2;
t89 = t87 * t81;
t90 = t86 * t82;
t70 = t80 * t89 - t90;
t99 = t57 * t70;
t98 = t57 * t70 ^ 2;
t96 = t61 * t67;
t75 = 0.1e1 / t79;
t77 = 0.1e1 / t81;
t95 = t75 * t77;
t93 = t79 * t87;
t71 = t80 * t92 + t91;
t66 = 0.1e1 / t71 ^ 2;
t88 = t66 * t100 * t87 ^ 2;
t78 = 0.1e1 / t81 ^ 2;
t76 = 0.1e1 / t100;
t69 = t80 * t90 - t89;
t65 = 0.1e1 / t71;
t64 = 0.1e1 / (0.1e1 + t88);
t62 = 0.1e1 / (t67 ^ 2 * t76 * t78 + 0.1e1);
t59 = t70 * t66 * t64 * t93;
t56 = 0.1e1 / t58;
t55 = (t67 * t76 * t77 * t80 + t86) * t62;
t54 = 0.1e1 / (0.1e1 + t98);
t53 = (t67 * t78 * t82 - t69 * t77) * t75 * t62;
t52 = (t71 * t56 - (t61 * t79 * t82 - t60 * t69 + (-t60 * t94 - t96) * t53) * t99) * t54;
t1 = [-t70 * t62 * t95, t55, 0, t53, t53, 0; (-t67 * t56 - (-t60 + (t95 * t96 + t60) * t62) * t98) * t54 (t55 * t96 * t99 + (-t56 * t93 - (t61 * t80 + (-t55 + t86) * t79 * t60) * t99) * t81) * t54, 0, t52, t52, 0; (-t66 * t69 * t87 + t65 * t86) * t79 * t64 (-t65 * t80 * t87 - t82 * t88) * t64, 0, -t59, -t59, 0;];
Ja_rot  = t1;
