% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PPRRRR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobia_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:42
% EndTime: 2019-02-26 19:42:43
% DurationCPUTime: 0.17s
% Computational Cost: add. (429->30), mult. (1141->67), div. (40->9), fcn. (1556->15), ass. (0->44)
t82 = sin(pkin(13));
t83 = sin(pkin(12));
t86 = cos(pkin(13));
t87 = cos(pkin(12));
t88 = cos(pkin(7));
t89 = cos(pkin(6));
t97 = t87 * t89;
t84 = sin(pkin(7));
t85 = sin(pkin(6));
t98 = t85 * t84;
t101 = (-t83 * t82 + t86 * t97) * t88 - t87 * t98;
t100 = t83 * t89;
t99 = t84 * t89;
t91 = cos(qJ(3));
t96 = t88 * t91;
t95 = t83 * t98;
t76 = -t86 * t100 - t87 * t82;
t77 = -t82 * t100 + t87 * t86;
t90 = sin(qJ(3));
t68 = t77 * t91 + (t76 * t88 + t95) * t90;
t72 = t83 * t85 * t88 - t76 * t84;
t81 = qJ(4) + qJ(5);
t79 = sin(t81);
t80 = cos(t81);
t59 = t68 * t80 + t72 * t79;
t57 = 0.1e1 / t59 ^ 2;
t58 = t68 * t79 - t72 * t80;
t93 = t58 ^ 2 * t57 + 0.1e1;
t75 = t82 * t97 + t83 * t86;
t71 = t90 * t99 + (t86 * t88 * t90 + t82 * t91) * t85;
t70 = -t91 * t99 + (t82 * t90 - t86 * t96) * t85;
t69 = 0.1e1 / t70 ^ 2;
t67 = -t76 * t96 + t77 * t90 - t91 * t95;
t66 = t101 * t90 + t75 * t91;
t64 = -t101 * t91 + t75 * t90;
t63 = atan2(-t64, t70);
t61 = cos(t63);
t60 = sin(t63);
t56 = 0.1e1 / t93;
t55 = -t60 * t64 + t61 * t70;
t54 = 0.1e1 / t55 ^ 2;
t52 = (-t66 / t70 + t71 * t64 * t69) / (t64 ^ 2 * t69 + 0.1e1);
t51 = t93 * t56;
t1 = [0, 0, t52, 0, 0, 0; 0, 0 (t68 / t55 - (-t60 * t66 + t61 * t71 + (-t60 * t70 - t61 * t64) * t52) * t67 * t54) / (t67 ^ 2 * t54 + 0.1e1) 0, 0, 0; 0, 0 (-t79 / t59 + t80 * t58 * t57) * t67 * t56, t51, t51, 0;];
Ja_rot  = t1;
