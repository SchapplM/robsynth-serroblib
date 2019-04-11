% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR10V2_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobia_rot_5_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:31
% EndTime: 2019-04-11 14:56:32
% DurationCPUTime: 0.16s
% Computational Cost: add. (600->33), mult. (873->84), div. (144->11), fcn. (1297->11), ass. (0->49)
t85 = qJ(2) + qJ(3);
t81 = sin(t85);
t87 = sin(qJ(4));
t100 = t81 * t87;
t82 = cos(t85);
t90 = cos(qJ(4));
t91 = cos(qJ(1));
t94 = t90 * t91;
t88 = sin(qJ(1));
t96 = t88 * t87;
t72 = t82 * t96 + t94;
t71 = atan2(-t72, t100);
t68 = sin(t71);
t69 = cos(t71);
t62 = t69 * t100 - t68 * t72;
t61 = 0.1e1 / t62 ^ 2;
t93 = t91 * t87;
t95 = t88 * t90;
t75 = t82 * t93 - t95;
t105 = t61 * t75;
t104 = t61 * t75 ^ 2;
t76 = t82 * t94 + t96;
t86 = sin(qJ(5));
t89 = cos(qJ(5));
t97 = t81 * t91;
t67 = t76 * t89 + t86 * t97;
t65 = 0.1e1 / t67 ^ 2;
t66 = t76 * t86 - t89 * t97;
t103 = t65 * t66;
t102 = t69 * t72;
t79 = 0.1e1 / t81;
t83 = 0.1e1 / t87;
t101 = t79 * t83;
t99 = t81 * t88;
t98 = t81 * t90;
t92 = t65 * t66 ^ 2 + 0.1e1;
t84 = 0.1e1 / t87 ^ 2;
t80 = 0.1e1 / t81 ^ 2;
t74 = t82 * t95 - t93;
t70 = 0.1e1 / (t72 ^ 2 * t80 * t84 + 0.1e1);
t64 = 0.1e1 / t67;
t63 = 0.1e1 / t92;
t60 = 0.1e1 / t62;
t59 = (t72 * t80 * t82 * t83 + t88) * t70;
t58 = 0.1e1 / (0.1e1 + t104);
t57 = (t72 * t84 * t90 - t74 * t83) * t79 * t70;
t56 = ((-t82 * t89 - t86 * t98) * t64 - (t82 * t86 - t89 * t98) * t103) * t63 * t91;
t55 = (t59 * t102 * t105 + (-t60 * t97 - (t69 * t82 + (-t59 * t81 + t99) * t68) * t105) * t87) * t58;
t1 = [-t75 * t70 * t101, t59, t59, t57, 0, 0; (-t72 * t60 - (-t68 + (t101 * t102 + t68) * t70) * t104) * t58, t55, t55 (t76 * t60 - (t69 * t98 - t68 * t74 + (-t68 * t100 - t102) * t57) * t105) * t58, 0, 0; ((-t74 * t86 + t89 * t99) * t64 - (-t74 * t89 - t86 * t99) * t103) * t63, t56, t56 (t89 * t103 - t64 * t86) * t75 * t63, t92 * t63, 0;];
Ja_rot  = t1;
