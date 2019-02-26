% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPPRR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:50
% EndTime: 2019-02-26 19:45:50
% DurationCPUTime: 0.22s
% Computational Cost: add. (632->36), mult. (1727->92), div. (65->9), fcn. (2425->15), ass. (0->56)
t89 = sin(pkin(6));
t97 = cos(qJ(5));
t104 = t89 * t97;
t92 = cos(pkin(6));
t98 = cos(qJ(2));
t102 = t92 * t98;
t88 = sin(pkin(10));
t91 = cos(pkin(10));
t95 = sin(qJ(2));
t82 = -t91 * t102 + t88 * t95;
t103 = t92 * t95;
t83 = t91 * t103 + t88 * t98;
t87 = sin(pkin(11));
t90 = cos(pkin(11));
t70 = t82 * t87 + t83 * t90;
t94 = sin(qJ(5));
t64 = -t91 * t104 + t70 * t94;
t81 = (-t87 * t98 + t90 * t95) * t89;
t77 = t81 * t94 + t92 * t97;
t63 = atan2(-t64, t77);
t60 = sin(t63);
t61 = cos(t63);
t54 = -t60 * t64 + t61 * t77;
t53 = 0.1e1 / t54 ^ 2;
t84 = t88 * t102 + t91 * t95;
t85 = -t88 * t103 + t91 * t98;
t74 = t84 * t87 + t85 * t90;
t67 = t88 * t104 + t74 * t94;
t109 = t53 * t67;
t100 = -t84 * t90 + t85 * t87;
t105 = t89 * t94;
t68 = -t88 * t105 + t74 * t97;
t93 = sin(qJ(6));
t96 = cos(qJ(6));
t59 = t100 * t93 + t68 * t96;
t57 = 0.1e1 / t59 ^ 2;
t58 = -t100 * t96 + t68 * t93;
t108 = t57 * t58;
t76 = 0.1e1 / t77 ^ 2;
t107 = t64 * t76;
t106 = t100 * t97;
t101 = t58 ^ 2 * t57 + 0.1e1;
t99 = -t60 * t77 - t61 * t64;
t80 = (t87 * t95 + t90 * t98) * t89;
t78 = t81 * t97 - t92 * t94;
t75 = 0.1e1 / t77;
t69 = -t82 * t90 + t83 * t87;
t66 = t91 * t105 + t70 * t97;
t62 = 0.1e1 / (t64 ^ 2 * t76 + 0.1e1);
t56 = 0.1e1 / t59;
t55 = 0.1e1 / t101;
t52 = 0.1e1 / t54;
t51 = 0.1e1 / (t67 ^ 2 * t53 + 0.1e1);
t50 = (t80 * t107 - t69 * t75) * t94 * t62;
t49 = (t78 * t107 - t66 * t75) * t62;
t1 = [0, t50, 0, 0, t49, 0; 0 (t100 * t94 * t52 - ((-t60 * t69 + t61 * t80) * t94 + t99 * t50) * t109) * t51, 0, 0 (t68 * t52 - (t99 * t49 - t60 * t66 + t61 * t78) * t109) * t51, 0; 0 ((t93 * t106 + t74 * t96) * t56 - (t96 * t106 - t74 * t93) * t108) * t55, 0, 0 (t96 * t108 - t56 * t93) * t67 * t55, t101 * t55;];
Ja_rot  = t1;
