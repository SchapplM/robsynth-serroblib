% Calculate inertial parameters regressor of potential energy for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR14V3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14V3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_energypot_fixb_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:10:05
% EndTime: 2019-04-12 15:10:05
% DurationCPUTime: 0.09s
% Computational Cost: add. (80->44), mult. (211->64), div. (0->0), fcn. (240->10), ass. (0->32)
t72 = sin(qJ(2));
t77 = cos(qJ(2));
t73 = sin(qJ(1));
t78 = cos(qJ(1));
t81 = g(1) * t78 + g(2) * t73;
t58 = -g(3) * t77 + t81 * t72;
t71 = sin(qJ(4));
t88 = t71 * t72;
t87 = t72 * t73;
t76 = cos(qJ(4));
t86 = t72 * t76;
t85 = t72 * t78;
t84 = t73 * t77;
t83 = t78 * t71;
t82 = t78 * t76;
t63 = t76 * t84 - t83;
t65 = t73 * t71 + t77 * t82;
t70 = sin(qJ(5));
t75 = cos(qJ(5));
t80 = g(1) * (t65 * t70 - t75 * t85) + g(2) * (t63 * t70 - t75 * t87) + g(3) * (t70 * t86 + t77 * t75);
t62 = t71 * t84 + t82;
t64 = -t73 * t76 + t77 * t83;
t79 = g(1) * t64 + g(2) * t62 + g(3) * t88;
t74 = cos(qJ(6));
t69 = sin(qJ(6));
t66 = g(1) * t73 - g(2) * t78;
t61 = -t77 * t70 + t75 * t86;
t59 = -g(3) * t72 - t81 * t77;
t57 = t58 * qJ(3);
t56 = t65 * t75 + t70 * t85;
t54 = t63 * t75 + t70 * t87;
t1 = [0, 0, 0, 0, 0, 0, -t81, t66, -g(3), 0, 0, 0, 0, 0, 0, 0, t59, t58, -t66, 0, 0, 0, 0, 0, 0, 0, t59, -t66, -t58, -t57, 0, 0, 0, 0, 0, 0, -g(1) * t65 - g(2) * t63 - g(3) * t86, t79, -t58, -t57, 0, 0, 0, 0, 0, 0, -g(1) * t56 - g(2) * t54 - g(3) * t61, t80, -t79, -t57, 0, 0, 0, 0, 0, 0, -g(1) * (t56 * t74 + t64 * t69) - g(2) * (t54 * t74 + t62 * t69) - g(3) * (t61 * t74 + t69 * t88) -g(1) * (-t56 * t69 + t64 * t74) - g(2) * (-t54 * t69 + t62 * t74) - g(3) * (-t61 * t69 + t74 * t88) -t80, -t57;];
U_reg  = t1;
