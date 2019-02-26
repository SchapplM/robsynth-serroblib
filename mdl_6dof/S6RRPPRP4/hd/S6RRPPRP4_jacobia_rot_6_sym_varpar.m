% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRP4_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:47
% EndTime: 2019-02-26 21:26:47
% DurationCPUTime: 0.21s
% Computational Cost: add. (456->29), mult. (1327->81), div. (77->9), fcn. (1878->11), ass. (0->49)
t89 = sin(qJ(1));
t91 = cos(qJ(2));
t100 = t89 * t91;
t86 = cos(pkin(9));
t85 = sin(pkin(9));
t92 = cos(qJ(1));
t99 = t92 * t85;
t82 = t100 * t86 - t99;
t87 = sin(qJ(5));
t90 = cos(qJ(5));
t98 = t92 * t86;
t93 = t100 * t85 + t98;
t71 = t82 * t90 + t87 * t93;
t110 = -t82 * t87 + t90 * t93;
t88 = sin(qJ(2));
t109 = t88 ^ 2;
t96 = t85 * t90 - t86 * t87;
t79 = t96 * t88;
t65 = atan2(t110, -t79);
t62 = sin(t65);
t63 = cos(t65);
t61 = t110 * t62 - t63 * t79;
t60 = 0.1e1 / t61 ^ 2;
t83 = t89 * t85 + t91 * t98;
t94 = -t89 * t86 + t91 * t99;
t72 = t83 * t87 - t90 * t94;
t108 = t60 * t72;
t107 = t63 * t110;
t73 = t83 * t90 + t87 * t94;
t68 = 0.1e1 / t73 ^ 2;
t106 = t68 * t92;
t78 = 0.1e1 / t79 ^ 2;
t105 = t110 * t78;
t104 = t72 ^ 2 * t60;
t101 = t88 * t92;
t97 = t62 * t79 + t107;
t95 = t85 * t87 + t86 * t90;
t81 = t96 * t91;
t80 = t95 * t88;
t77 = 0.1e1 / t79;
t74 = t89 * t79;
t67 = 0.1e1 / t73;
t66 = 0.1e1 / (t109 * t68 * t92 ^ 2 + 0.1e1);
t64 = 0.1e1 / (t110 ^ 2 * t78 + 0.1e1);
t59 = 0.1e1 / t61;
t58 = 0.1e1 / (0.1e1 + t104);
t57 = (t105 * t81 + t74 * t77) * t64;
t56 = (-t105 * t80 + t71 * t77) * t64;
t1 = [t72 * t77 * t64, t57, 0, 0, t56, 0; (t110 * t59 - (-t62 + (t107 * t77 + t62) * t64) * t104) * t58 (-(t57 * t97 - t62 * t74 - t63 * t81) * t108 + t96 * t59 * t101) * t58, 0, 0 (t73 * t59 - (t56 * t97 - t62 * t71 + t63 * t80) * t108) * t58, 0; (t71 * t106 - t89 * t67) * t88 * t66 (t106 * t109 * t95 + t91 * t67) * t66 * t92, 0, 0, t72 * t68 * t66 * t101, 0;];
Ja_rot  = t1;
