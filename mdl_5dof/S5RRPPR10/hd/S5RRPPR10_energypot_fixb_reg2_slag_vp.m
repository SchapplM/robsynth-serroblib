% Calculate inertial parameters regressor of potential energy for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:39
% EndTime: 2019-12-31 19:44:40
% DurationCPUTime: 0.12s
% Computational Cost: add. (111->59), mult. (231->79), div. (0->0), fcn. (250->8), ass. (0->36)
t75 = cos(pkin(8));
t77 = sin(qJ(2));
t95 = t75 * t77;
t74 = sin(pkin(8));
t96 = t74 * t77;
t101 = pkin(3) * t95 + qJ(4) * t96;
t100 = g(3) * pkin(5);
t99 = pkin(7) * t77;
t98 = g(3) * t77;
t97 = t77 * pkin(2) + pkin(5);
t78 = sin(qJ(1));
t80 = cos(qJ(2));
t94 = t78 * t80;
t81 = cos(qJ(1));
t93 = t81 * t74;
t92 = t81 * t75;
t91 = t81 * pkin(1) + t78 * pkin(6);
t90 = qJ(3) * t77;
t89 = t78 * pkin(1) - t81 * pkin(6);
t88 = t91 + (pkin(2) * t80 + t90) * t81;
t87 = -t80 * qJ(3) + t97;
t86 = g(1) * t81 + g(2) * t78;
t85 = pkin(2) * t94 + t78 * t90 + t89;
t59 = -t78 * t75 + t80 * t93;
t60 = t78 * t74 + t80 * t92;
t84 = t60 * pkin(3) + t59 * qJ(4) + t88;
t57 = t74 * t94 + t92;
t83 = g(1) * t59 + g(2) * t57 + g(3) * t96;
t58 = t75 * t94 - t93;
t82 = t58 * pkin(3) + t57 * qJ(4) + t85;
t79 = cos(qJ(5));
t76 = sin(qJ(5));
t61 = g(1) * t78 - g(2) * t81;
t54 = -g(3) * t80 + t86 * t77;
t53 = -g(1) * t60 - g(2) * t58 - g(3) * t95;
t1 = [0, 0, 0, 0, 0, 0, -t86, t61, -g(3), -t100, 0, 0, 0, 0, 0, 0, -t86 * t80 - t98, t54, -t61, -g(1) * t91 - g(2) * t89 - t100, 0, 0, 0, 0, 0, 0, t53, t83, -t54, -g(1) * t88 - g(2) * t85 - g(3) * t87, 0, 0, 0, 0, 0, 0, t53, -t54, -t83, -g(1) * t84 - g(2) * t82 - g(3) * (t87 + t101), 0, 0, 0, 0, 0, 0, -g(1) * (t59 * t76 + t60 * t79) - g(2) * (t57 * t76 + t58 * t79) - (t74 * t76 + t75 * t79) * t98, -g(1) * (t59 * t79 - t60 * t76) - g(2) * (t57 * t79 - t58 * t76) - (t74 * t79 - t75 * t76) * t98, t54, -g(1) * (t60 * pkin(4) - t81 * t99 + t84) - g(2) * (t58 * pkin(4) - t78 * t99 + t82) - g(3) * (pkin(4) * t95 + (pkin(7) - qJ(3)) * t80 + t97 + t101);];
U_reg = t1;
