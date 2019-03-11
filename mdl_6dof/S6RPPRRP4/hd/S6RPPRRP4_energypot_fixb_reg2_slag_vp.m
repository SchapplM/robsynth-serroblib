% Calculate inertial parameters regressor of potential energy for
% S6RPPRRP4
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
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:06:18
% EndTime: 2019-03-09 02:06:18
% DurationCPUTime: 0.14s
% Computational Cost: add. (158->55), mult. (286->69), div. (0->0), fcn. (340->8), ass. (0->38)
t79 = sin(qJ(4));
t81 = cos(qJ(4));
t105 = -pkin(4) * t81 - pkin(8) * t79;
t104 = g(3) * pkin(6);
t77 = -qJ(3) + pkin(6);
t101 = g(3) * t77;
t100 = g(3) * t79;
t99 = sin(qJ(1));
t78 = sin(qJ(5));
t98 = t78 * t81;
t80 = cos(qJ(5));
t97 = t80 * t81;
t82 = cos(qJ(1));
t96 = t82 * pkin(1) + t99 * qJ(2);
t95 = cos(pkin(9));
t94 = sin(pkin(9));
t93 = t81 * pkin(8) + t77;
t92 = t82 * pkin(2) + t96;
t91 = t99 * pkin(1) - t82 * qJ(2);
t90 = t99 * pkin(2) + t91;
t64 = -t82 * t95 - t99 * t94;
t65 = t82 * t94 - t99 * t95;
t89 = g(1) * t65 - g(2) * t64;
t88 = g(1) * t64 + g(2) * t65;
t87 = -t64 * pkin(3) + t65 * pkin(7) + t92;
t86 = -t65 * pkin(3) - t64 * pkin(7) + t90;
t52 = t64 * t80 - t65 * t98;
t54 = -t64 * t98 - t65 * t80;
t85 = -g(1) * t54 - g(2) * t52 + t78 * t100;
t84 = t105 * t64 + t87;
t83 = t105 * t65 + t86;
t67 = -g(1) * t82 - g(2) * t99;
t66 = g(1) * t99 - g(2) * t82;
t55 = -t64 * t97 + t65 * t78;
t53 = -t64 * t78 - t65 * t97;
t51 = g(3) * t81 - t88 * t79;
t50 = -g(1) * t55 - g(2) * t53 + t80 * t100;
t1 = [0, 0, 0, 0, 0, 0, t67, t66, -g(3), -t104, 0, 0, 0, 0, 0, 0, t67, -g(3), -t66, -g(1) * t96 - g(2) * t91 - t104, 0, 0, 0, 0, 0, 0, t88, t89, g(3), -g(1) * t92 - g(2) * t90 - t101, 0, 0, 0, 0, 0, 0, t88 * t81 + t100, t51, -t89, -g(1) * t87 - g(2) * t86 - t101, 0, 0, 0, 0, 0, 0, t50, -t85, -t51, -g(1) * t84 - g(2) * t83 - g(3) * (-t79 * pkin(4) + t93) 0, 0, 0, 0, 0, 0, t50, -t51, t85, -g(1) * (t55 * pkin(5) + t54 * qJ(6) + t84) - g(2) * (t53 * pkin(5) + t52 * qJ(6) + t83) - g(3) * t93 - (-pkin(5) * t80 - qJ(6) * t78 - pkin(4)) * t100;];
U_reg  = t1;
