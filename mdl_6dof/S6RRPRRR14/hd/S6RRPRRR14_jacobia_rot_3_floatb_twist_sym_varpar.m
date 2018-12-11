% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR14
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR14_jacobia_rot_3_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_rot_3_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_rot_3_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:20
% EndTime: 2018-12-10 18:38:21
% DurationCPUTime: 0.21s
% Computational Cost: add. (825->43), mult. (979->85), div. (50->9), fcn. (1135->21), ass. (0->58)
t93 = sin(pkin(6));
t96 = sin(qJ(1));
t106 = t93 * t96;
t90 = pkin(6) + qJ(2);
t84 = cos(t90) / 0.2e1;
t103 = pkin(6) - qJ(2);
t87 = cos(t103);
t78 = t87 / 0.2e1 + t84;
t95 = sin(qJ(2));
t98 = cos(qJ(1));
t69 = -t96 * t78 - t98 * t95;
t100 = sin(t103);
t102 = sin(t90) / 0.2e1;
t76 = t102 - t100 / 0.2e1;
t97 = cos(qJ(2));
t71 = t96 * t76 - t98 * t97;
t88 = pkin(7) + pkin(14);
t80 = sin(t88) / 0.2e1;
t89 = pkin(7) - pkin(14);
t85 = sin(t89);
t73 = t80 - t85 / 0.2e1;
t81 = cos(t89) / 0.2e1;
t86 = cos(t88);
t75 = t81 - t86 / 0.2e1;
t94 = cos(pkin(14));
t55 = t75 * t106 + t69 * t73 - t71 * t94;
t53 = 0.1e1 / t55 ^ 2;
t72 = t80 + t85 / 0.2e1;
t74 = t81 + t86 / 0.2e1;
t91 = sin(pkin(14));
t54 = -t72 * t106 - t69 * t74 - t71 * t91;
t109 = t53 * t54;
t104 = cos(pkin(7));
t101 = t93 * t104;
t66 = -t98 * t78 + t96 * t95;
t92 = sin(pkin(7));
t60 = t98 * t101 - t66 * t92;
t65 = -(t102 + t100 / 0.2e1) * t92 + cos(pkin(6)) * t104;
t59 = atan2(t60, t65);
t57 = cos(t59);
t108 = t57 * t60;
t56 = sin(t59);
t50 = t56 * t60 + t57 * t65;
t49 = 0.1e1 / t50 ^ 2;
t61 = -t96 * t101 + t69 * t92;
t107 = t61 ^ 2 * t49;
t105 = t93 * t98;
t99 = -t98 * t76 - t96 * t97;
t77 = t84 - t87 / 0.2e1;
t64 = 0.1e1 / t65 ^ 2;
t63 = 0.1e1 / t65;
t58 = 0.1e1 / (t60 ^ 2 * t64 + 0.1e1);
t52 = 0.1e1 / t55;
t51 = 0.1e1 / (t54 ^ 2 * t53 + 0.1e1);
t48 = 0.1e1 / t50;
t47 = 0.1e1 / (0.1e1 + t107);
t46 = (t60 * t64 * t77 + t63 * t99) * t92 * t58;
t1 = [t61 * t63 * t58, t46, 0, 0, 0, 0; (t60 * t48 + (t56 + (t63 * t108 - t56) * t58) * t107) * t47 (-t71 * t92 * t48 + ((t56 * t99 - t57 * t77) * t92 + (-t56 * t65 + t108) * t46) * t61 * t49) * t47, 0, 0, 0, 0; ((-t72 * t105 - t66 * t74 + t91 * t99) * t52 - (t75 * t105 + t66 * t73 + t94 * t99) * t109) * t51 ((t69 * t91 - t71 * t74) * t52 - (t69 * t94 + t71 * t73) * t109) * t51, 0, 0, 0, 0;];
Ja_rot  = t1;
