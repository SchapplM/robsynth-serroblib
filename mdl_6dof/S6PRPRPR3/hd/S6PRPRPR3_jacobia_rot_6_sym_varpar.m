% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:38
% EndTime: 2019-02-26 19:47:39
% DurationCPUTime: 0.21s
% Computational Cost: add. (632->32), mult. (1727->83), div. (65->9), fcn. (2425->15), ass. (0->55)
t87 = sin(pkin(6));
t92 = sin(qJ(4));
t103 = t87 * t92;
t90 = cos(pkin(6));
t85 = sin(pkin(11));
t88 = cos(pkin(11));
t93 = sin(qJ(2));
t96 = cos(qJ(2));
t98 = t96 * t85 + t93 * t88;
t80 = t98 * t90;
t81 = t93 * t85 - t96 * t88;
t86 = sin(pkin(10));
t89 = cos(pkin(10));
t69 = t89 * t80 - t86 * t81;
t95 = cos(qJ(4));
t64 = -t89 * t103 + t69 * t95;
t79 = t98 * t87;
t76 = t79 * t95 + t90 * t92;
t62 = atan2(-t64, t76);
t59 = sin(t62);
t60 = cos(t62);
t53 = -t59 * t64 + t60 * t76;
t52 = 0.1e1 / t53 ^ 2;
t99 = -t86 * t80 - t89 * t81;
t67 = t86 * t103 + t95 * t99;
t108 = t52 * t67;
t97 = t81 * t90;
t71 = t86 * t97 - t89 * t98;
t94 = cos(qJ(6));
t104 = t71 * t94;
t102 = t87 * t95;
t66 = -t86 * t102 + t92 * t99;
t91 = sin(qJ(6));
t58 = t66 * t91 - t104;
t56 = 0.1e1 / t58 ^ 2;
t105 = t71 * t91;
t57 = -t66 * t94 - t105;
t107 = t56 * t57;
t74 = 0.1e1 / t76 ^ 2;
t106 = t64 * t74;
t101 = t57 ^ 2 * t56 + 0.1e1;
t100 = -t59 * t76 - t60 * t64;
t78 = t81 * t87;
t75 = -t79 * t92 + t90 * t95;
t73 = 0.1e1 / t76;
t68 = -t86 * t98 - t89 * t97;
t63 = t89 * t102 + t69 * t92;
t61 = 0.1e1 / (t64 ^ 2 * t74 + 0.1e1);
t55 = 0.1e1 / t58;
t54 = 0.1e1 / t101;
t51 = 0.1e1 / t53;
t50 = 0.1e1 / (t67 ^ 2 * t52 + 0.1e1);
t49 = (-t78 * t106 - t68 * t73) * t95 * t61;
t48 = (t75 * t106 + t63 * t73) * t61;
t1 = [0, t49, 0, t48, 0, 0; 0 (t71 * t95 * t51 - ((-t59 * t68 - t60 * t78) * t95 + t100 * t49) * t108) * t50, 0 (-t66 * t51 - (t100 * t48 + t59 * t63 + t60 * t75) * t108) * t50, 0, 0; 0 ((-t92 * t104 + t91 * t99) * t55 - (t92 * t105 + t94 * t99) * t107) * t54, 0 (-t91 * t107 - t55 * t94) * t67 * t54, 0, t101 * t54;];
Ja_rot  = t1;
