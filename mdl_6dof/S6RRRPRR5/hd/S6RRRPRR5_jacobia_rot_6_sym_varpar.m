% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR5_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:17
% EndTime: 2019-02-26 22:18:17
% DurationCPUTime: 0.10s
% Computational Cost: add. (419->22), mult. (359->54), div. (84->9), fcn. (524->9), ass. (0->40)
t94 = qJ(2) + qJ(3);
t92 = cos(t94);
t95 = sin(qJ(1));
t100 = t95 * t92;
t90 = sin(t94);
t84 = atan2(-t100, t90);
t82 = sin(t84);
t83 = cos(t84);
t75 = -t82 * t100 + t83 * t90;
t74 = 0.1e1 / t75 ^ 2;
t96 = cos(qJ(1));
t108 = t74 * t96 ^ 2;
t93 = qJ(5) + qJ(6);
t91 = cos(t93);
t101 = t95 * t91;
t89 = sin(t93);
t99 = t96 * t89;
t81 = t90 * t99 + t101;
t79 = 0.1e1 / t81 ^ 2;
t102 = t95 * t89;
t98 = t96 * t91;
t80 = -t90 * t98 + t102;
t107 = t79 * t80;
t106 = t82 * t90;
t88 = t92 ^ 2;
t105 = 0.1e1 / t90 ^ 2 * t88;
t104 = t92 * t96;
t85 = 0.1e1 / (t95 ^ 2 * t105 + 0.1e1);
t103 = t95 * t85;
t97 = t80 ^ 2 * t79 + 0.1e1;
t86 = 0.1e1 / t90;
t78 = 0.1e1 / t81;
t77 = (0.1e1 + t105) * t103;
t76 = 0.1e1 / t97;
t73 = 0.1e1 / t75;
t72 = 0.1e1 / (t88 * t108 + 0.1e1);
t71 = t97 * t76;
t70 = (-t89 * t107 - t78 * t91) * t76 * t104;
t69 = (-t90 * t73 - (t95 * t106 + t83 * t92 + (-t83 * t100 - t106) * t77) * t92 * t74) * t96 * t72;
t1 = [-t86 * t85 * t104, t77, t77, 0, 0, 0; (-t73 * t100 - (t83 * t86 * t88 * t103 + (t85 - 0.1e1) * t92 * t82) * t92 * t108) * t72, t69, t69, 0, 0, 0; ((t90 * t101 + t99) * t78 - (-t90 * t102 + t98) * t107) * t76, t70, t70, 0, t71, t71;];
Ja_rot  = t1;
