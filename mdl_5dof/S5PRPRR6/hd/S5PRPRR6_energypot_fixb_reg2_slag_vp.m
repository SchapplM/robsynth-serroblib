% Calculate inertial parameters regressor of potential energy for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:52
% EndTime: 2019-12-05 15:57:52
% DurationCPUTime: 0.22s
% Computational Cost: add. (207->81), mult. (360->122), div. (0->0), fcn. (422->12), ass. (0->44)
t80 = sin(pkin(9));
t81 = sin(pkin(5));
t105 = t80 * t81;
t83 = cos(pkin(9));
t104 = t81 * t83;
t87 = sin(qJ(2));
t103 = t81 * t87;
t89 = cos(qJ(2));
t102 = t81 * t89;
t79 = sin(pkin(10));
t84 = cos(pkin(5));
t101 = t84 * t79;
t100 = t84 * t87;
t99 = t84 * t89;
t98 = t83 * pkin(1) + pkin(6) * t105;
t97 = t84 * pkin(6) + qJ(1);
t96 = t79 * t105;
t75 = t80 * pkin(1);
t95 = -pkin(6) * t104 + t75;
t62 = t80 * t99 + t83 * t87;
t63 = -t80 * t100 + t83 * t89;
t82 = cos(pkin(10));
t72 = t82 * pkin(3) + pkin(2);
t85 = -pkin(7) - qJ(3);
t94 = pkin(3) * t96 - t62 * t85 + t63 * t72 + t98;
t93 = g(1) * t80 - g(2) * t83;
t92 = pkin(3) * t101 + t85 * t102 + t72 * t103 + t97;
t61 = t83 * t100 + t80 * t89;
t78 = pkin(10) + qJ(4);
t73 = sin(t78);
t74 = cos(t78);
t50 = t74 * t104 + t61 * t73;
t52 = -t74 * t105 + t63 * t73;
t56 = t73 * t103 - t84 * t74;
t91 = g(1) * t52 + g(2) * t50 + g(3) * t56;
t60 = t80 * t87 - t83 * t99;
t49 = -g(1) * t62 - g(2) * t60 + g(3) * t102;
t90 = t61 * t72 - t60 * t85 + t75 + (-pkin(3) * t79 - pkin(6)) * t104;
t88 = cos(qJ(5));
t86 = sin(qJ(5));
t57 = t74 * t103 + t84 * t73;
t53 = t73 * t105 + t63 * t74;
t51 = -t73 * t104 + t61 * t74;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t80, t93, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t63 - g(2) * t61 - g(3) * t103, -t49, -g(3) * t84 - t93 * t81, -g(1) * t98 - g(2) * t95 - g(3) * t97, 0, 0, 0, 0, 0, 0, -g(1) * (t63 * t82 + t96) - g(2) * (-t79 * t104 + t61 * t82) - g(3) * (t82 * t103 + t101), -g(1) * (t82 * t105 - t63 * t79) - g(2) * (-t82 * t104 - t61 * t79) - g(3) * (-t79 * t103 + t84 * t82), t49, -g(1) * (t63 * pkin(2) + t62 * qJ(3) + t98) - g(2) * (t61 * pkin(2) + t60 * qJ(3) + t95) - g(3) * ((pkin(2) * t87 - qJ(3) * t89) * t81 + t97), 0, 0, 0, 0, 0, 0, -g(1) * t53 - g(2) * t51 - g(3) * t57, t91, t49, -g(1) * t94 - g(2) * t90 - g(3) * t92, 0, 0, 0, 0, 0, 0, -g(1) * (t53 * t88 + t62 * t86) - g(2) * (t51 * t88 + t60 * t86) - g(3) * (-t86 * t102 + t57 * t88), -g(1) * (-t53 * t86 + t62 * t88) - g(2) * (-t51 * t86 + t60 * t88) - g(3) * (-t88 * t102 - t57 * t86), -t91, -g(1) * (t53 * pkin(4) + t52 * pkin(8) + t94) - g(2) * (t51 * pkin(4) + t50 * pkin(8) + t90) - g(3) * (t57 * pkin(4) + t56 * pkin(8) + t92);];
U_reg = t1;
