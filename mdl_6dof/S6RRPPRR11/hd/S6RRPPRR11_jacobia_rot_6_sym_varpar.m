% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR11_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:19
% EndTime: 2019-02-26 21:34:19
% DurationCPUTime: 0.16s
% Computational Cost: add. (878->37), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
t80 = cos(pkin(6));
t85 = cos(qJ(2));
t86 = cos(qJ(1));
t91 = t86 * t85;
t82 = sin(qJ(2));
t83 = sin(qJ(1));
t94 = t83 * t82;
t72 = -t80 * t91 + t94;
t78 = pkin(11) + qJ(5);
t76 = sin(t78);
t77 = cos(t78);
t79 = sin(pkin(6));
t95 = t79 * t86;
t64 = t72 * t77 + t76 * t95;
t96 = t79 * t85;
t70 = t80 * t76 + t77 * t96;
t61 = atan2(t64, t70);
t54 = sin(t61);
t55 = cos(t61);
t52 = t54 * t64 + t55 * t70;
t51 = 0.1e1 / t52 ^ 2;
t92 = t86 * t82;
t93 = t83 * t85;
t87 = t80 * t93 + t92;
t97 = t79 * t83;
t62 = t76 * t97 - t87 * t77;
t105 = t51 * t62;
t104 = t55 * t64;
t74 = -t80 * t94 + t91;
t81 = sin(qJ(6));
t100 = t74 * t81;
t63 = t87 * t76 + t77 * t97;
t84 = cos(qJ(6));
t60 = t63 * t84 + t100;
t57 = 0.1e1 / t60 ^ 2;
t99 = t74 * t84;
t59 = t63 * t81 - t99;
t103 = t57 * t59;
t102 = t62 ^ 2 * t51;
t69 = 0.1e1 / t70 ^ 2;
t101 = t64 * t69;
t98 = t79 * t82;
t90 = t59 ^ 2 * t57 + 0.1e1;
t89 = -t54 * t70 + t104;
t88 = -t72 * t76 + t77 * t95;
t73 = t80 * t92 + t93;
t71 = -t76 * t96 + t80 * t77;
t68 = 0.1e1 / t70;
t58 = 0.1e1 / (t64 ^ 2 * t69 + 0.1e1);
t56 = 0.1e1 / t60;
t53 = 0.1e1 / t90;
t50 = 0.1e1 / t52;
t49 = 0.1e1 / (0.1e1 + t102);
t48 = (t98 * t101 + t68 * t73) * t77 * t58;
t47 = (-t71 * t101 + t68 * t88) * t58;
t1 = [-t62 * t68 * t58, t48, 0, 0, t47, 0; (t64 * t50 - (-t54 + (-t68 * t104 + t54) * t58) * t102) * t49 (-t74 * t77 * t50 - ((t54 * t73 - t55 * t98) * t77 + t89 * t48) * t105) * t49, 0, 0 (t63 * t50 - (t89 * t47 + t54 * t88 + t55 * t71) * t105) * t49, 0; ((t73 * t84 + t81 * t88) * t56 - (-t73 * t81 + t84 * t88) * t103) * t53 ((t76 * t100 + t87 * t84) * t56 - (t76 * t99 - t87 * t81) * t103) * t53, 0, 0 (t84 * t103 - t81 * t56) * t62 * t53, t90 * t53;];
Ja_rot  = t1;
