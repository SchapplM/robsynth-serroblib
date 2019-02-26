% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR7_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobia_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:29
% EndTime: 2019-02-26 19:57:30
% DurationCPUTime: 0.17s
% Computational Cost: add. (563->40), mult. (1633->92), div. (35->9), fcn. (2228->17), ass. (0->52)
t100 = cos(qJ(2));
t93 = cos(pkin(13));
t105 = t93 * t100;
t96 = cos(pkin(6));
t98 = sin(qJ(2));
t107 = t96 * t98;
t88 = sin(pkin(13));
t85 = -t88 * t107 + t105;
t87 = sin(pkin(14));
t112 = t85 * t87;
t90 = sin(pkin(7));
t91 = sin(pkin(6));
t111 = t90 * t91;
t94 = cos(pkin(8));
t110 = t90 * t94;
t95 = cos(pkin(7));
t109 = t91 * t95;
t92 = cos(pkin(14));
t108 = t92 * t95;
t106 = t88 * t100;
t84 = -t96 * t106 - t93 * t98;
t101 = t88 * t111 + t84 * t95;
t75 = t101 * t92 - t112;
t81 = t88 * t109 - t84 * t90;
t89 = sin(pkin(8));
t103 = t75 * t94 + t81 * t89;
t76 = t101 * t87 + t85 * t92;
t97 = sin(qJ(4));
t99 = cos(qJ(4));
t65 = t103 * t97 + t76 * t99;
t63 = 0.1e1 / t65 ^ 2;
t64 = -t103 * t99 + t76 * t97;
t104 = t64 ^ 2 * t63 + 0.1e1;
t77 = -t85 * t108 - t84 * t87;
t102 = t85 * t89 * t90 + t77 * t94;
t83 = t93 * t107 + t106;
t82 = t96 * t105 - t88 * t98;
t79 = (-(-t100 * t87 - t98 * t108) * t89 + t98 * t110) * t91;
t78 = -t95 * t112 + t84 * t92;
t74 = -(t96 * t90 * t92 + (t100 * t108 - t87 * t98) * t91) * t89 + (-t100 * t111 + t96 * t95) * t94;
t73 = 0.1e1 / t74 ^ 2;
t72 = (-t83 * t108 - t82 * t87) * t89 - t83 * t110;
t71 = -t75 * t89 + t81 * t94;
t70 = (-t83 * t87 + (-t93 * t111 + t82 * t95) * t92) * t89 - (-t93 * t109 - t82 * t90) * t94;
t69 = atan2(t70, t74);
t67 = cos(t69);
t66 = sin(t69);
t62 = 0.1e1 / t104;
t61 = t66 * t70 + t67 * t74;
t60 = 0.1e1 / t61 ^ 2;
t58 = (t72 / t74 - t79 * t70 * t73) / (t70 ^ 2 * t73 + 0.1e1);
t1 = [0, t58, 0, 0, 0, 0; 0 ((t85 * t110 - t77 * t89) / t61 - (t66 * t72 + t67 * t79 + (-t66 * t74 + t67 * t70) * t58) * t71 * t60) / (t71 ^ 2 * t60 + 0.1e1) 0, 0, 0, 0; 0 ((-t102 * t99 + t78 * t97) / t65 - (t102 * t97 + t78 * t99) * t64 * t63) * t62, 0, t104 * t62, 0, 0;];
Ja_rot  = t1;
