% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_jacobia_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:26:20
% EndTime: 2019-02-26 22:26:20
% DurationCPUTime: 0.11s
% Computational Cost: add. (469->27), mult. (656->69), div. (146->11), fcn. (1013->9), ass. (0->42)
t84 = qJ(2) + qJ(3);
t79 = cos(t84);
t85 = sin(qJ(4));
t88 = cos(qJ(1));
t90 = t88 * t85;
t86 = sin(qJ(1));
t87 = cos(qJ(4));
t91 = t86 * t87;
t71 = t79 * t91 - t90;
t78 = sin(t84);
t93 = t78 * t87;
t69 = atan2(-t71, t93);
t65 = sin(t69);
t66 = cos(t69);
t64 = -t65 * t71 + t66 * t93;
t63 = 0.1e1 / t64 ^ 2;
t89 = t88 * t87;
t92 = t86 * t85;
t74 = t79 * t89 + t92;
t100 = t63 * t74;
t98 = t66 * t71;
t73 = -t79 * t90 + t91;
t77 = 0.1e1 / t78 ^ 2;
t83 = 0.1e1 / t88 ^ 2;
t68 = 0.1e1 / (t73 ^ 2 * t83 * t77 + 0.1e1);
t76 = 0.1e1 / t78;
t97 = t68 * t76;
t96 = t74 ^ 2 * t63;
t80 = 0.1e1 / t87;
t95 = t76 * t80;
t94 = t77 * t79;
t82 = 0.1e1 / t88;
t81 = 0.1e1 / t87 ^ 2;
t70 = t79 * t92 + t89;
t67 = 0.1e1 / (t71 ^ 2 * t77 * t81 + 0.1e1);
t62 = 0.1e1 / t64;
t61 = (t71 * t80 * t94 + t86) * t67;
t60 = (-t73 * t82 * t94 + t85) * t68;
t59 = 0.1e1 / (0.1e1 + t96);
t58 = (-t71 * t81 * t85 + t70 * t80) * t76 * t67;
t57 = (t61 * t98 * t100 + (-t88 * t78 * t62 - (t66 * t79 + (-t61 + t86) * t78 * t65) * t100) * t87) * t59;
t1 = [-t74 * t67 * t95, t61, t61, t58, 0, 0; (-t71 * t62 - (-t65 + (t95 * t98 + t65) * t67) * t96) * t59, t57, t57 (t73 * t62 - (-t66 * t78 * t85 + t65 * t70 + (-t65 * t93 - t98) * t58) * t100) * t59, 0, 0; (t73 * t83 * t86 + t70 * t82) * t97, t60, t60, -t74 * t82 * t97, 0, 0;];
Ja_rot  = t1;
