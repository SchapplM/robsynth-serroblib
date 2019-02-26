% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR4_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:39:32
% EndTime: 2019-02-26 21:39:32
% DurationCPUTime: 0.13s
% Computational Cost: add. (422->27), mult. (1051->69), div. (55->9), fcn. (1490->13), ass. (0->46)
t83 = sin(qJ(1));
t85 = cos(qJ(1));
t79 = sin(pkin(11));
t81 = cos(pkin(11));
t84 = cos(qJ(2));
t92 = cos(pkin(6));
t89 = t84 * t92;
t82 = sin(qJ(2));
t90 = t82 * t92;
t86 = -t79 * t90 + t81 * t89;
t87 = t84 * t79 + t82 * t81;
t62 = -t83 * t87 + t85 * t86;
t73 = t82 * t79 - t84 * t81;
t80 = sin(pkin(6));
t70 = t73 * t80;
t56 = atan2(t62, t70);
t54 = cos(t56);
t97 = t54 * t62;
t72 = t79 * t89 + t81 * t90;
t66 = -t83 * t72 - t85 * t73;
t78 = qJ(4) + pkin(12);
t76 = sin(t78);
t77 = cos(t78);
t94 = t80 * t83;
t60 = t66 * t77 + t76 * t94;
t58 = 0.1e1 / t60 ^ 2;
t59 = t66 * t76 - t77 * t94;
t96 = t58 * t59;
t53 = sin(t56);
t51 = t53 * t62 + t54 * t70;
t50 = 0.1e1 / t51 ^ 2;
t64 = -t83 * t86 - t85 * t87;
t95 = t64 ^ 2 * t50;
t93 = t80 * t85;
t91 = t59 ^ 2 * t58 + 0.1e1;
t88 = -t85 * t72 + t83 * t73;
t71 = t87 * t80;
t69 = 0.1e1 / t70 ^ 2;
t68 = 0.1e1 / t70;
t57 = 0.1e1 / t60;
t55 = 0.1e1 / (t62 ^ 2 * t69 + 0.1e1);
t52 = 0.1e1 / t91;
t49 = 0.1e1 / t51;
t48 = 0.1e1 / (0.1e1 + t95);
t47 = (-t62 * t69 * t71 + t68 * t88) * t55;
t1 = [t64 * t68 * t55, t47, 0, 0, 0, 0; (t62 * t49 + (t53 + (t68 * t97 - t53) * t55) * t95) * t48 (t66 * t49 + (t53 * t88 + t54 * t71 + (-t53 * t70 + t97) * t47) * t64 * t50) * t48, 0, 0, 0, 0; ((t76 * t88 - t77 * t93) * t57 - (t76 * t93 + t77 * t88) * t96) * t52 (t76 * t57 - t77 * t96) * t64 * t52, 0, t91 * t52, 0, 0;];
Ja_rot  = t1;
