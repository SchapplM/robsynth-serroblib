% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR6_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:56:58
% EndTime: 2019-02-26 21:56:58
% DurationCPUTime: 0.22s
% Computational Cost: add. (994->25), mult. (1174->55), div. (115->9), fcn. (1719->11), ass. (0->38)
t103 = qJ(4) + qJ(5);
t100 = cos(t103);
t111 = sin(qJ(2));
t89 = sin(t103);
t93 = cos(qJ(2));
t83 = t111 * t100 - t93 * t89;
t91 = sin(qJ(1));
t75 = t83 * t91;
t82 = t93 * t100 + t111 * t89;
t70 = atan2(t75, t82);
t68 = cos(t70);
t108 = t68 * t75;
t81 = 0.1e1 / t82 ^ 2;
t69 = 0.1e1 / (t75 ^ 2 * t81 + 0.1e1);
t77 = t82 * t91;
t80 = 0.1e1 / t82;
t113 = (-t75 * t81 * t83 - t77 * t80) * t69;
t67 = sin(t70);
t65 = t67 * t75 + t68 * t82;
t64 = 0.1e1 / t65 ^ 2;
t112 = cos(qJ(1));
t79 = t83 * t112;
t105 = t79 ^ 2 * t64;
t62 = 0.1e1 / (0.1e1 + t105);
t63 = 0.1e1 / t65;
t78 = t82 * t112;
t58 = ((t113 * (-t67 * t82 + t108) - t67 * t77 + t68 * t83) * t64 * t79 + t78 * t63) * t62;
t90 = sin(qJ(6));
t92 = cos(qJ(6));
t74 = t78 * t92 - t91 * t90;
t72 = 0.1e1 / t74 ^ 2;
t73 = t78 * t90 + t91 * t92;
t107 = t72 * t73;
t102 = t73 ^ 2 * t72 + 0.1e1;
t66 = 0.1e1 / t102;
t71 = 0.1e1 / t74;
t59 = (t92 * t107 - t90 * t71) * t66 * t79;
t1 = [t79 * t80 * t69, -t113, 0, t113, t113, 0; (t75 * t63 - (-t67 + (-t80 * t108 + t67) * t69) * t105) * t62, -t58, 0, t58, t58, 0; ((t112 * t92 - t77 * t90) * t71 - (-t112 * t90 - t77 * t92) * t107) * t66, t59, 0, -t59, -t59, t102 * t66;];
Ja_rot  = t1;
