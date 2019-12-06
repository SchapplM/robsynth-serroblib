% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:40
% EndTime: 2019-12-05 16:19:40
% DurationCPUTime: 0.12s
% Computational Cost: add. (107->30), mult. (277->70), div. (0->0), fcn. (184->6), ass. (0->28)
t114 = qJD(2) ^ 2;
t120 = t114 / 0.2e1;
t113 = cos(qJ(3));
t116 = qJD(2) * t113;
t111 = sin(qJ(3));
t118 = pkin(6) * t116 + t111 * qJD(1);
t100 = qJ(4) * t116 + t118;
t108 = sin(pkin(9));
t109 = cos(pkin(9));
t106 = t113 * qJD(1);
t117 = qJD(2) * t111;
t98 = qJD(3) * pkin(3) + t106 + (-pkin(6) - qJ(4)) * t117;
t91 = t109 * t100 + t108 * t98;
t119 = t113 * t114;
t115 = qJD(2) * qJD(3);
t90 = -t100 * t108 + t109 * t98;
t103 = qJD(4) + (-pkin(3) * t113 - pkin(2)) * qJD(2);
t112 = cos(qJ(5));
t110 = sin(qJ(5));
t107 = qJD(3) + qJD(5);
t102 = (t108 * t113 + t109 * t111) * qJD(2);
t101 = (-t108 * t111 + t109 * t113) * qJD(2);
t94 = -t101 * pkin(4) + t103;
t93 = t101 * t110 + t102 * t112;
t92 = -t112 * t101 + t102 * t110;
t89 = pkin(7) * t101 + t91;
t88 = qJD(3) * pkin(4) - pkin(7) * t102 + t90;
t1 = [qJD(1) ^ 2 / 0.2e1, t120, 0, 0, t111 ^ 2 * t120, t111 * t119, t111 * t115, t113 * t115, qJD(3) ^ 2 / 0.2e1, pkin(2) * t119 + (-pkin(6) * t117 + t106) * qJD(3), -t114 * pkin(2) * t111 - qJD(3) * t118, t101 * t91 - t102 * t90, t91 ^ 2 / 0.2e1 + t90 ^ 2 / 0.2e1 + t103 ^ 2 / 0.2e1, t93 ^ 2 / 0.2e1, -t93 * t92, t93 * t107, -t92 * t107, t107 ^ 2 / 0.2e1, t94 * t92 + (-t110 * t89 + t112 * t88) * t107, t94 * t93 - (t110 * t88 + t112 * t89) * t107;];
T_reg = t1;
