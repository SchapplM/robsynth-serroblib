% Calculate inertial parameters regressor of potential energy for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:09:36
% EndTime: 2019-03-09 03:09:36
% DurationCPUTime: 0.13s
% Computational Cost: add. (227->69), mult. (201->88), div. (0->0), fcn. (201->10), ass. (0->41)
t87 = qJ(2) + pkin(6);
t110 = g(3) * t87;
t89 = sin(qJ(3));
t109 = g(3) * t89;
t84 = qJ(1) + pkin(9);
t77 = sin(t84);
t85 = sin(pkin(10));
t108 = t77 * t85;
t91 = cos(qJ(3));
t107 = t77 * t91;
t79 = cos(t84);
t106 = t79 * t91;
t105 = t85 * t91;
t86 = cos(pkin(10));
t104 = t86 * t91;
t88 = -pkin(8) - qJ(4);
t103 = t88 * t89;
t90 = sin(qJ(1));
t102 = t90 * pkin(1) + t77 * pkin(2);
t92 = cos(qJ(1));
t101 = t92 * pkin(1) + t79 * pkin(2) + t77 * pkin(7);
t74 = t86 * pkin(4) + pkin(3);
t100 = t89 * t74 + t91 * t88 + t87;
t99 = -t79 * pkin(7) + t102;
t98 = g(1) * t79 + g(2) * t77;
t97 = -g(1) * t92 - g(2) * t90;
t96 = pkin(3) * t91 + qJ(4) * t89;
t83 = pkin(10) + qJ(5);
t76 = sin(t83);
t78 = cos(t83);
t60 = t76 * t107 + t79 * t78;
t62 = t76 * t106 - t77 * t78;
t95 = g(1) * t62 + g(2) * t60 + t76 * t109;
t94 = pkin(4) * t108 - t79 * t103 + t74 * t106 + t101;
t93 = -t77 * t103 + t74 * t107 + (-pkin(4) * t85 - pkin(7)) * t79 + t102;
t67 = g(1) * t77 - g(2) * t79;
t64 = -g(3) * t91 + t98 * t89;
t63 = t78 * t106 + t77 * t76;
t61 = t78 * t107 - t79 * t76;
t59 = -g(1) * t63 - g(2) * t61 - t78 * t109;
t1 = [0, 0, 0, 0, 0, 0, t97, g(1) * t90 - g(2) * t92, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t98, t67, -g(3), t97 * pkin(1) - t110, 0, 0, 0, 0, 0, 0, -t98 * t91 - t109, t64, -t67, -g(1) * t101 - g(2) * t99 - t110, 0, 0, 0, 0, 0, 0, -g(1) * (t79 * t104 + t108) - g(2) * (t77 * t104 - t79 * t85) - t86 * t109, -g(1) * (-t79 * t105 + t77 * t86) - g(2) * (-t77 * t105 - t79 * t86) + t85 * t109, -t64, -g(1) * (t96 * t79 + t101) - g(2) * (t96 * t77 + t99) - g(3) * (t89 * pkin(3) - t91 * qJ(4) + t87) 0, 0, 0, 0, 0, 0, t59, t95, -t64, -g(1) * t94 - g(2) * t93 - g(3) * t100, 0, 0, 0, 0, 0, 0, t59, -t64, -t95, -g(1) * (t63 * pkin(5) + t62 * qJ(6) + t94) - g(2) * (t61 * pkin(5) + t60 * qJ(6) + t93) - g(3) * ((pkin(5) * t78 + qJ(6) * t76) * t89 + t100);];
U_reg  = t1;
