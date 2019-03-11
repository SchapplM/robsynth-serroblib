% Calculate inertial parameters regressor of potential energy for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:13:00
% EndTime: 2019-03-09 03:13:00
% DurationCPUTime: 0.14s
% Computational Cost: add. (186->61), mult. (188->69), div. (0->0), fcn. (184->8), ass. (0->39)
t86 = sin(qJ(3));
t101 = qJ(4) * t86;
t83 = qJ(1) + pkin(9);
t77 = sin(t83);
t89 = cos(qJ(3));
t106 = t77 * t89;
t109 = pkin(3) * t106 + t77 * t101;
t84 = qJ(2) + pkin(6);
t108 = g(3) * t84;
t107 = g(3) * t89;
t78 = cos(t83);
t105 = t78 * t89;
t85 = sin(qJ(5));
t104 = t85 * t86;
t88 = cos(qJ(5));
t103 = t86 * t88;
t87 = sin(qJ(1));
t102 = t87 * pkin(1) + t77 * pkin(2);
t100 = t86 * pkin(3) + t84;
t90 = cos(qJ(1));
t99 = t90 * pkin(1) + t78 * pkin(2) + t77 * pkin(7);
t98 = -t78 * pkin(7) + t102;
t97 = pkin(3) * t105 + t78 * t101 + t99;
t96 = g(1) * t78 + g(2) * t77;
t95 = -g(1) * t90 - g(2) * t87;
t94 = -t89 * qJ(4) + t100;
t93 = t77 * pkin(4) + pkin(8) * t105 + t97;
t60 = -t78 * t103 + t77 * t85;
t62 = t77 * t103 + t78 * t85;
t92 = g(1) * t60 - g(2) * t62 + t88 * t107;
t91 = pkin(8) * t106 + (-pkin(4) - pkin(7)) * t78 + t102 + t109;
t79 = t86 * pkin(8);
t64 = g(1) * t77 - g(2) * t78;
t63 = t77 * t104 - t78 * t88;
t61 = t78 * t104 + t77 * t88;
t59 = g(3) * t86 + t96 * t89;
t58 = t96 * t86 - t107;
t57 = -g(1) * t61 - g(2) * t63 + t85 * t107;
t1 = [0, 0, 0, 0, 0, 0, t95, g(1) * t87 - g(2) * t90, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t96, t64, -g(3), t95 * pkin(1) - t108, 0, 0, 0, 0, 0, 0, -t59, t58, -t64, -g(1) * t99 - g(2) * t98 - t108, 0, 0, 0, 0, 0, 0, -t64, t59, -t58, -g(1) * t97 - g(2) * (t98 + t109) - g(3) * t94, 0, 0, 0, 0, 0, 0, t57, t92, -t59, -g(1) * t93 - g(2) * t91 - g(3) * (t79 + t94) 0, 0, 0, 0, 0, 0, t57, -t59, -t92, -g(1) * (t61 * pkin(5) + t60 * qJ(6) + t93) - g(2) * (t63 * pkin(5) - t62 * qJ(6) + t91) - g(3) * (t79 + t100) - (-pkin(5) * t85 + qJ(6) * t88 - qJ(4)) * t107;];
U_reg  = t1;
