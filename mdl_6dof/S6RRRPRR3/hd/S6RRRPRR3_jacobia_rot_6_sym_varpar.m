% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR3
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
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR3_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:07
% EndTime: 2019-02-26 22:17:07
% DurationCPUTime: 0.21s
% Computational Cost: add. (994->26), mult. (1174->55), div. (115->9), fcn. (1719->11), ass. (0->38)
t102 = sin(qJ(1));
t101 = sin(qJ(5));
t113 = qJ(2) + qJ(3);
t110 = sin(t113);
t121 = cos(qJ(5));
t99 = cos(t113);
t93 = -t99 * t101 + t110 * t121;
t85 = t93 * t102;
t92 = t110 * t101 + t99 * t121;
t80 = atan2(t85, t92);
t78 = cos(t80);
t118 = t78 * t85;
t91 = 0.1e1 / t92 ^ 2;
t79 = 0.1e1 / (t85 ^ 2 * t91 + 0.1e1);
t87 = t92 * t102;
t90 = 0.1e1 / t92;
t123 = (-t85 * t91 * t93 - t87 * t90) * t79;
t77 = sin(t80);
t75 = t77 * t85 + t78 * t92;
t74 = 0.1e1 / t75 ^ 2;
t122 = cos(qJ(1));
t89 = t93 * t122;
t115 = t89 ^ 2 * t74;
t72 = 0.1e1 / (0.1e1 + t115);
t73 = 0.1e1 / t75;
t88 = t92 * t122;
t68 = (((-t77 * t92 + t118) * t123 - t77 * t87 + t78 * t93) * t74 * t89 + t88 * t73) * t72;
t100 = sin(qJ(6));
t103 = cos(qJ(6));
t84 = -t102 * t100 + t88 * t103;
t82 = 0.1e1 / t84 ^ 2;
t83 = t88 * t100 + t102 * t103;
t117 = t82 * t83;
t112 = t83 ^ 2 * t82 + 0.1e1;
t76 = 0.1e1 / t112;
t81 = 0.1e1 / t84;
t69 = (-t100 * t81 + t103 * t117) * t76 * t89;
t1 = [t89 * t90 * t79, -t123, -t123, 0, t123, 0; (t85 * t73 - (-t77 + (-t90 * t118 + t77) * t79) * t115) * t72, -t68, -t68, 0, t68, 0; ((-t87 * t100 + t122 * t103) * t81 - (-t122 * t100 - t87 * t103) * t117) * t76, t69, t69, 0, -t69, t112 * t76;];
Ja_rot  = t1;
