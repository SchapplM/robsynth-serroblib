% Calculate minimal parameter regressor of potential energy for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:47:31
% EndTime: 2021-01-15 21:47:31
% DurationCPUTime: 0.07s
% Computational Cost: add. (64->33), mult. (76->49), div. (0->0), fcn. (81->10), ass. (0->26)
t96 = qJ(2) + pkin(9);
t92 = sin(t96);
t114 = g(3) * t92;
t101 = sin(qJ(1));
t97 = qJ(4) + qJ(5);
t94 = sin(t97);
t113 = t101 * t94;
t95 = cos(t97);
t112 = t101 * t95;
t99 = sin(qJ(4));
t111 = t101 * t99;
t104 = cos(qJ(1));
t110 = t104 * t94;
t109 = t104 * t95;
t108 = t104 * t99;
t102 = cos(qJ(4));
t107 = t101 * t102;
t106 = t104 * t102;
t105 = g(1) * t104 + g(2) * t101;
t103 = cos(qJ(2));
t100 = sin(qJ(2));
t98 = pkin(6) + qJ(3);
t93 = cos(t96);
t91 = t103 * pkin(2) + pkin(1);
t90 = g(1) * t101 - g(2) * t104;
t1 = [0, -t105, t90, 0, 0, 0, 0, 0, -g(3) * t100 - t105 * t103, -g(3) * t103 + t105 * t100, -t105 * t93 - t114, -g(3) * t93 + t105 * t92, -t90, -g(1) * (t98 * t101 + t104 * t91) - g(2) * (t101 * t91 - t104 * t98) - g(3) * (t100 * pkin(2) + pkin(5)), 0, 0, 0, 0, 0, -g(1) * (t93 * t106 + t111) - g(2) * (t93 * t107 - t108) - t102 * t114, -g(1) * (-t93 * t108 + t107) - g(2) * (-t93 * t111 - t106) + t99 * t114, 0, 0, 0, 0, 0, -g(1) * (t93 * t109 + t113) - g(2) * (t93 * t112 - t110) - t95 * t114, -g(1) * (-t93 * t110 + t112) - g(2) * (-t93 * t113 - t109) + t94 * t114;];
U_reg = t1;
