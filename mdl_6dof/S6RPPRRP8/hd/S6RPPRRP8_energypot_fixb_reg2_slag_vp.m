% Calculate inertial parameters regressor of potential energy for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:16:23
% EndTime: 2019-03-09 02:16:23
% DurationCPUTime: 0.15s
% Computational Cost: add. (148->61), mult. (181->69), div. (0->0), fcn. (177->8), ass. (0->41)
t100 = g(3) * pkin(6);
t99 = pkin(2) + pkin(6);
t73 = sin(pkin(9));
t98 = pkin(3) * t73;
t72 = pkin(9) + qJ(4);
t66 = sin(t72);
t97 = pkin(4) * t66;
t67 = cos(t72);
t96 = pkin(8) * t67;
t95 = g(3) * t67;
t76 = sin(qJ(5));
t77 = sin(qJ(1));
t94 = t77 * t76;
t78 = cos(qJ(5));
t93 = t77 * t78;
t79 = cos(qJ(1));
t92 = t79 * t76;
t91 = t79 * t78;
t90 = t79 * pkin(1) + t77 * qJ(2);
t74 = cos(pkin(9));
t89 = t74 * pkin(3) + t99;
t88 = t77 * t98 + t90;
t87 = -qJ(2) - t98;
t70 = t77 * pkin(1);
t75 = -pkin(7) - qJ(3);
t86 = -t77 * t75 + t70;
t85 = -t79 * qJ(2) + t70;
t84 = t67 * pkin(4) + t66 * pkin(8) + t89;
t83 = t79 * t96 + t86;
t58 = g(1) * t77 - g(2) * t79;
t82 = t88 + (-t96 + t97) * t77;
t54 = t66 * t94 - t91;
t56 = t66 * t92 + t93;
t81 = g(1) * t54 - g(2) * t56 + t76 * t95;
t80 = (g(1) * t75 - g(2) * (t87 - t97)) * t79;
t59 = g(1) * t79 + g(2) * t77;
t57 = -t66 * t91 + t94;
t55 = t66 * t93 + t92;
t53 = -g(3) * t66 + t58 * t67;
t52 = -g(1) * t55 - g(2) * t57 - t78 * t95;
t1 = [0, 0, 0, 0, 0, 0, -t59, t58, -g(3), -t100, 0, 0, 0, 0, 0, 0, -g(3), t59, -t58, -g(1) * t90 - g(2) * t85 - t100, 0, 0, 0, 0, 0, 0, -g(3) * t74 - t58 * t73, g(3) * t73 - t58 * t74, -t59, -g(1) * (t79 * qJ(3) + t90) - g(2) * (t77 * qJ(3) + t85) - g(3) * t99, 0, 0, 0, 0, 0, 0, -t58 * t66 - t95, -t53, -t59, -g(1) * (-t79 * t75 + t88) - g(2) * (t87 * t79 + t86) - g(3) * t89, 0, 0, 0, 0, 0, 0, t52, t81, t53, -g(1) * t82 - g(2) * t83 - g(3) * t84 + t80, 0, 0, 0, 0, 0, 0, t52, t53, -t81, -g(1) * (t55 * pkin(5) + t54 * qJ(6) + t82) - g(2) * (t57 * pkin(5) - t56 * qJ(6) + t83) - g(3) * ((pkin(5) * t78 + qJ(6) * t76) * t67 + t84) + t80;];
U_reg  = t1;
