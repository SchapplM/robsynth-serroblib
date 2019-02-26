% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR11_jacobia_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobia_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobia_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:31
% EndTime: 2019-02-26 20:54:31
% DurationCPUTime: 0.19s
% Computational Cost: add. (507->39), mult. (1523->86), div. (50->9), fcn. (2097->15), ass. (0->56)
t84 = sin(pkin(7));
t88 = cos(pkin(7));
t83 = sin(pkin(12));
t92 = cos(qJ(1));
t101 = t92 * t83;
t89 = cos(pkin(6));
t110 = sin(qJ(1));
t87 = cos(pkin(12));
t96 = t110 * t87;
t94 = t89 * t96 + t101;
t85 = sin(pkin(6));
t98 = t85 * t110;
t111 = -t84 * t98 + t94 * t88;
t100 = t92 * t87;
t97 = t110 * t83;
t79 = -t89 * t97 + t100;
t90 = sin(qJ(3));
t91 = cos(qJ(3));
t68 = -t111 * t90 + t79 * t91;
t74 = t84 * t94 + t88 * t98;
t82 = sin(pkin(13));
t86 = cos(pkin(13));
t58 = t68 * t86 + t74 * t82;
t56 = 0.1e1 / t58 ^ 2;
t57 = t68 * t82 - t74 * t86;
t109 = t56 * t57;
t102 = t88 * t91;
t78 = t101 * t89 + t96;
t106 = t78 * t90;
t77 = -t100 * t89 + t97;
t104 = t85 * t92;
t99 = t84 * t104;
t63 = t102 * t77 + t91 * t99 + t106;
t105 = t84 * t89;
t71 = -t91 * t105 + (-t102 * t87 + t83 * t90) * t85;
t62 = atan2(-t63, t71);
t60 = cos(t62);
t108 = t60 * t63;
t59 = sin(t62);
t53 = -t59 * t63 + t60 * t71;
t52 = 0.1e1 / t53 ^ 2;
t67 = t111 * t91 + t79 * t90;
t107 = t67 ^ 2 * t52;
t103 = t88 * t90;
t66 = t103 * t77 - t78 * t91 + t90 * t99;
t73 = t104 * t88 - t77 * t84;
t72 = t90 * t105 + (t103 * t87 + t83 * t91) * t85;
t70 = 0.1e1 / t71 ^ 2;
t69 = 0.1e1 / t71;
t61 = 0.1e1 / (t63 ^ 2 * t70 + 0.1e1);
t55 = 0.1e1 / t58;
t54 = 0.1e1 / (t57 ^ 2 * t56 + 0.1e1);
t51 = 0.1e1 / t53;
t50 = 0.1e1 / (0.1e1 + t107);
t49 = (t63 * t70 * t72 + t66 * t69) * t61;
t1 = [-t67 * t69 * t61, 0, t49, 0, 0, 0; ((-t106 + (-t77 * t88 - t99) * t91) * t51 - (-t59 + (t108 * t69 + t59) * t61) * t107) * t50, 0 (t68 * t51 - (t59 * t66 + t60 * t72 + (-t59 * t71 - t108) * t49) * t67 * t52) * t50, 0, 0, 0; ((t66 * t82 - t73 * t86) * t55 - (t66 * t86 + t73 * t82) * t109) * t54, 0 (t109 * t86 - t82 * t55) * t67 * t54, 0, 0, 0;];
Ja_rot  = t1;
