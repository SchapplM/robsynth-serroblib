% Calculate minimal parameter regressor of potential energy for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:55
% EndTime: 2019-12-05 17:39:55
% DurationCPUTime: 0.04s
% Computational Cost: add. (46->22), mult. (56->28), div. (0->0), fcn. (50->8), ass. (0->15)
t128 = sin(qJ(1));
t129 = cos(qJ(1));
t131 = t129 * pkin(1) + t128 * qJ(2);
t125 = pkin(8) + qJ(4);
t130 = t128 * pkin(1) - t129 * qJ(2);
t115 = g(1) * t128 - g(2) * t129;
t127 = cos(pkin(8));
t126 = sin(pkin(8));
t121 = qJ(5) + t125;
t120 = cos(t125);
t119 = sin(t125);
t118 = cos(t121);
t117 = sin(t121);
t116 = g(1) * t129 + g(2) * t128;
t1 = [0, -t116, t115, t116, -t115, -g(3) * pkin(5) - g(1) * t131 - g(2) * t130, -g(3) * t127 - t115 * t126, g(3) * t126 - t115 * t127, -t116, -g(1) * (t129 * qJ(3) + t131) - g(2) * (t128 * qJ(3) + t130) - g(3) * (pkin(2) + pkin(5)), 0, 0, 0, 0, 0, -g(3) * t120 - t115 * t119, g(3) * t119 - t115 * t120, 0, 0, 0, 0, 0, -g(3) * t118 - t115 * t117, g(3) * t117 - t115 * t118;];
U_reg = t1;
