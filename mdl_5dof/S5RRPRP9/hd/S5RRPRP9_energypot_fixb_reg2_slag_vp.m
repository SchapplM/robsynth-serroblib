% Calculate inertial parameters regressor of potential energy for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:27
% EndTime: 2019-12-31 20:07:27
% DurationCPUTime: 0.14s
% Computational Cost: add. (129->61), mult. (184->77), div. (0->0), fcn. (187->8), ass. (0->37)
t98 = g(3) * pkin(5);
t78 = sin(qJ(2));
t97 = g(3) * t78;
t77 = -pkin(7) - qJ(3);
t96 = t77 * t78;
t75 = sin(pkin(8));
t79 = sin(qJ(1));
t95 = t79 * t75;
t80 = cos(qJ(2));
t94 = t79 * t80;
t74 = pkin(8) + qJ(4);
t68 = sin(t74);
t81 = cos(qJ(1));
t93 = t81 * t68;
t69 = cos(t74);
t92 = t81 * t69;
t91 = t81 * t75;
t76 = cos(pkin(8));
t90 = t81 * t76;
t89 = t81 * pkin(1) + t79 * pkin(6);
t66 = t76 * pkin(3) + pkin(2);
t88 = t78 * t66 + t80 * t77 + pkin(5);
t71 = t79 * pkin(1);
t87 = -t81 * pkin(6) + t71;
t86 = g(1) * t81 + g(2) * t79;
t85 = pkin(2) * t80 + qJ(3) * t78;
t84 = pkin(3) * t95 + t89 + (t66 * t80 - t96) * t81;
t55 = t68 * t94 + t92;
t57 = -t79 * t69 + t80 * t93;
t83 = g(1) * t57 + g(2) * t55 + t68 * t97;
t82 = -t79 * t96 + t66 * t94 + t71 + (-pkin(3) * t75 - pkin(6)) * t81;
t63 = g(1) * t79 - g(2) * t81;
t59 = -g(3) * t80 + t78 * t86;
t58 = t79 * t68 + t80 * t92;
t56 = t69 * t94 - t93;
t54 = -g(1) * t58 - g(2) * t56 - t69 * t97;
t1 = [0, 0, 0, 0, 0, 0, -t86, t63, -g(3), -t98, 0, 0, 0, 0, 0, 0, -t80 * t86 - t97, t59, -t63, -g(1) * t89 - g(2) * t87 - t98, 0, 0, 0, 0, 0, 0, -g(1) * (t80 * t90 + t95) - g(2) * (t76 * t94 - t91) - t76 * t97, -g(1) * (t79 * t76 - t80 * t91) - g(2) * (-t75 * t94 - t90) + t75 * t97, -t59, -g(1) * (t81 * t85 + t89) - g(2) * (t85 * t79 + t87) - g(3) * (t78 * pkin(2) - t80 * qJ(3) + pkin(5)), 0, 0, 0, 0, 0, 0, t54, t83, -t59, -g(1) * t84 - g(2) * t82 - g(3) * t88, 0, 0, 0, 0, 0, 0, t54, -t59, -t83, -g(1) * (t58 * pkin(4) + t57 * qJ(5) + t84) - g(2) * (t56 * pkin(4) + t55 * qJ(5) + t82) - g(3) * ((pkin(4) * t69 + qJ(5) * t68) * t78 + t88);];
U_reg = t1;
