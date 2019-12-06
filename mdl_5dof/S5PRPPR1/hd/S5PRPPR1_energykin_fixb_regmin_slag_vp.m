% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:08
% EndTime: 2019-12-05 15:22:08
% DurationCPUTime: 0.13s
% Computational Cost: add. (100->31), mult. (254->73), div. (0->0), fcn. (157->6), ass. (0->24)
t106 = sin(pkin(8));
t108 = cos(pkin(8));
t112 = t108 * qJD(2);
t100 = qJ(3) * t112 + t106 * qJD(1);
t105 = sin(pkin(9));
t107 = cos(pkin(9));
t96 = qJD(3) + (-pkin(3) * t108 - qJ(4) * t106 - pkin(2)) * qJD(2);
t90 = t107 * t100 + t105 * t96;
t114 = t106 * t107;
t113 = qJD(2) * t106;
t111 = t105 * t113;
t89 = -t105 * t100 + t107 * t96;
t99 = -qJ(3) * t113 + t108 * qJD(1);
t98 = qJD(4) - t99;
t110 = cos(qJ(5));
t109 = sin(qJ(5));
t104 = -qJD(2) * pkin(2) + qJD(3);
t101 = -qJD(5) + t112;
t93 = (-t105 * t109 + t107 * t110) * t113;
t92 = (t105 * t110 + t107 * t109) * t113;
t91 = pkin(4) * t111 + t98;
t88 = -pkin(6) * t111 + t90;
t87 = (-pkin(4) * t108 - pkin(6) * t114) * qJD(2) + t89;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, -t104 * t112, t104 * t113, (t100 * t108 - t106 * t99) * qJD(2), t100 ^ 2 / 0.2e1 + t99 ^ 2 / 0.2e1 + t104 ^ 2 / 0.2e1, (t105 * t106 * t98 - t108 * t89) * qJD(2), (t108 * t90 + t98 * t114) * qJD(2), (-t105 * t90 - t107 * t89) * t113, t90 ^ 2 / 0.2e1 + t89 ^ 2 / 0.2e1 + t98 ^ 2 / 0.2e1, t93 ^ 2 / 0.2e1, -t93 * t92, -t93 * t101, t92 * t101, t101 ^ 2 / 0.2e1, -(-t109 * t88 + t110 * t87) * t101 + t91 * t92, (t109 * t87 + t110 * t88) * t101 + t91 * t93;];
T_reg = t1;
