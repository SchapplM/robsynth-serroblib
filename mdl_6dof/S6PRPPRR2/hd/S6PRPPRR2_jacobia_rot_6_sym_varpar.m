% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPPRR2_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobia_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:12
% EndTime: 2019-02-26 19:45:13
% DurationCPUTime: 0.22s
% Computational Cost: add. (632->31), mult. (1727->84), div. (65->9), fcn. (2425->15), ass. (0->54)
t89 = sin(pkin(6));
t94 = sin(qJ(5));
t104 = t89 * t94;
t87 = sin(pkin(11));
t90 = cos(pkin(11));
t95 = sin(qJ(2));
t98 = cos(qJ(2));
t85 = t95 * t87 - t98 * t90;
t92 = cos(pkin(6));
t83 = t85 * t92;
t88 = sin(pkin(10));
t91 = cos(pkin(10));
t99 = t98 * t87 + t95 * t90;
t73 = t91 * t83 + t88 * t99;
t97 = cos(qJ(5));
t69 = t91 * t104 + t73 * t97;
t81 = t85 * t89;
t79 = -t81 * t97 + t92 * t94;
t66 = atan2(t69, t79);
t63 = sin(t66);
t64 = cos(t66);
t57 = t63 * t69 + t64 * t79;
t56 = 0.1e1 / t57 ^ 2;
t75 = t88 * t83 - t91 * t99;
t67 = t88 * t104 + t75 * t97;
t108 = t56 * t67;
t84 = t99 * t92;
t100 = -t88 * t84 - t91 * t85;
t103 = t89 * t97;
t68 = t88 * t103 - t75 * t94;
t93 = sin(qJ(6));
t96 = cos(qJ(6));
t62 = t100 * t93 + t68 * t96;
t60 = 0.1e1 / t62 ^ 2;
t61 = -t100 * t96 + t68 * t93;
t107 = t60 * t61;
t78 = 0.1e1 / t79 ^ 2;
t106 = t69 * t78;
t105 = t100 * t94;
t102 = t61 ^ 2 * t60 + 0.1e1;
t101 = -t63 * t79 + t64 * t69;
t82 = t99 * t89;
t80 = t81 * t94 + t92 * t97;
t77 = 0.1e1 / t79;
t72 = t91 * t84 - t88 * t85;
t70 = t91 * t103 - t73 * t94;
t65 = 0.1e1 / (t69 ^ 2 * t78 + 0.1e1);
t59 = 0.1e1 / t62;
t58 = 0.1e1 / t102;
t55 = 0.1e1 / t57;
t54 = 0.1e1 / (t67 ^ 2 * t56 + 0.1e1);
t53 = (t82 * t106 + t72 * t77) * t97 * t65;
t52 = (-t80 * t106 + t70 * t77) * t65;
t1 = [0, t53, 0, 0, t52, 0; 0 (-t100 * t97 * t55 - ((t63 * t72 - t64 * t82) * t97 + t101 * t53) * t108) * t54, 0, 0 (t68 * t55 - (t101 * t52 + t63 * t70 + t64 * t80) * t108) * t54, 0; 0 ((t93 * t105 - t75 * t96) * t59 - (t96 * t105 + t75 * t93) * t107) * t58, 0, 0 (t96 * t107 - t59 * t93) * t67 * t58, t102 * t58;];
Ja_rot  = t1;
