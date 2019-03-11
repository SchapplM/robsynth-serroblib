% Calculate inertial parameters regressor of potential energy for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:58:42
% EndTime: 2019-03-09 01:58:42
% DurationCPUTime: 0.13s
% Computational Cost: add. (201->62), mult. (160->73), div. (0->0), fcn. (152->10), ass. (0->38)
t84 = cos(qJ(5));
t67 = t84 * pkin(5) + pkin(4);
t75 = pkin(10) + qJ(4);
t68 = sin(t75);
t70 = cos(t75);
t79 = -qJ(6) - pkin(8);
t101 = t67 * t70 - t68 * t79;
t100 = g(3) * t68;
t80 = qJ(2) + pkin(6);
t99 = g(3) * t80;
t76 = qJ(1) + pkin(9);
t69 = sin(t76);
t82 = sin(qJ(5));
t96 = t69 * t82;
t95 = t69 * t84;
t71 = cos(t76);
t94 = t71 * t82;
t93 = t71 * t84;
t78 = cos(pkin(10));
t66 = t78 * pkin(3) + pkin(2);
t85 = cos(qJ(1));
t74 = t85 * pkin(1);
t92 = t71 * t66 + t74;
t83 = sin(qJ(1));
t73 = t83 * pkin(1);
t81 = -pkin(7) - qJ(3);
t91 = t69 * t66 + t71 * t81 + t73;
t77 = sin(pkin(10));
t90 = t77 * pkin(3) + t80;
t89 = -t69 * t81 + t92;
t88 = pkin(4) * t70 + pkin(8) * t68;
t87 = g(1) * t71 + g(2) * t69;
t86 = -g(1) * t85 - g(2) * t83;
t60 = g(1) * t69 - g(2) * t71;
t59 = -g(3) * t70 + t87 * t68;
t58 = -g(1) * (t70 * t93 + t96) - g(2) * (t70 * t95 - t94) - t84 * t100;
t57 = -g(1) * (-t70 * t94 + t95) - g(2) * (-t70 * t96 - t93) + t82 * t100;
t1 = [0, 0, 0, 0, 0, 0, t86, g(1) * t83 - g(2) * t85, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t87, t60, -g(3), t86 * pkin(1) - t99, 0, 0, 0, 0, 0, 0, -g(3) * t77 - t87 * t78, -g(3) * t78 + t87 * t77, -t60, -g(1) * (t71 * pkin(2) + t69 * qJ(3) + t74) - g(2) * (t69 * pkin(2) - t71 * qJ(3) + t73) - t99, 0, 0, 0, 0, 0, 0, -t87 * t70 - t100, t59, -t60, -g(1) * t89 - g(2) * t91 - g(3) * t90, 0, 0, 0, 0, 0, 0, t58, t57, -t59, -g(1) * (t88 * t71 + t89) - g(2) * (t88 * t69 + t91) - g(3) * (t68 * pkin(4) - t70 * pkin(8) + t90) 0, 0, 0, 0, 0, 0, t58, t57, -t59, -g(1) * (t101 * t71 + t92) - g(2) * (-pkin(5) * t94 + t91) - g(3) * (t68 * t67 + t70 * t79 + t90) + (-g(1) * (pkin(5) * t82 - t81) - g(2) * t101) * t69;];
U_reg  = t1;
