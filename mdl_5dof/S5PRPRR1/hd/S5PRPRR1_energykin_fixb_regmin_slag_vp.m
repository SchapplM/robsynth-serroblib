% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:03
% EndTime: 2019-12-05 15:43:03
% DurationCPUTime: 0.06s
% Computational Cost: add. (102->33), mult. (273->69), div. (0->0), fcn. (190->6), ass. (0->26)
t115 = cos(qJ(4));
t109 = sin(qJ(4));
t107 = cos(pkin(9));
t104 = t107 * qJD(1);
t106 = sin(pkin(9));
t113 = qJD(2) * t106;
t93 = t104 + (-pkin(6) - qJ(3)) * t113;
t112 = qJD(2) * t107;
t98 = qJ(3) * t112 + t106 * qJD(1);
t94 = pkin(6) * t112 + t98;
t114 = t109 * t93 + t115 * t94;
t111 = -t109 * t94 + t115 * t93;
t99 = qJD(3) + (-pkin(3) * t107 - pkin(2)) * qJD(2);
t110 = cos(qJ(5));
t108 = sin(qJ(5));
t105 = qJD(4) + qJD(5);
t102 = -qJD(2) * pkin(2) + qJD(3);
t97 = -qJ(3) * t113 + t104;
t96 = (t115 * t106 + t107 * t109) * qJD(2);
t95 = t109 * t113 - t115 * t112;
t88 = t95 * pkin(4) + t99;
t87 = -t108 * t95 + t110 * t96;
t86 = t108 * t96 + t110 * t95;
t85 = -t95 * pkin(7) + t114;
t84 = qJD(4) * pkin(4) - t96 * pkin(7) + t111;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, -t102 * t112, t102 * t113, (-t106 * t97 + t107 * t98) * qJD(2), t98 ^ 2 / 0.2e1 + t97 ^ 2 / 0.2e1 + t102 ^ 2 / 0.2e1, t96 ^ 2 / 0.2e1, -t96 * t95, t96 * qJD(4), -t95 * qJD(4), qJD(4) ^ 2 / 0.2e1, t111 * qJD(4) + t99 * t95, -t114 * qJD(4) + t99 * t96, t87 ^ 2 / 0.2e1, -t87 * t86, t87 * t105, -t86 * t105, t105 ^ 2 / 0.2e1, t88 * t86 + (-t108 * t85 + t110 * t84) * t105, t88 * t87 - (t108 * t84 + t110 * t85) * t105;];
T_reg = t1;
