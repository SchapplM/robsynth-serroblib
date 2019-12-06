% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPPR2
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
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:41
% EndTime: 2019-12-05 15:24:41
% DurationCPUTime: 0.06s
% Computational Cost: add. (80->28), mult. (187->64), div. (0->0), fcn. (120->8), ass. (0->27)
t108 = sin(pkin(9));
t110 = cos(pkin(9));
t115 = cos(qJ(2));
t102 = qJD(2) * pkin(2) + t115 * qJD(1);
t109 = sin(pkin(8));
t111 = cos(pkin(8));
t113 = sin(qJ(2));
t120 = qJD(1) * t113;
t98 = t109 * t102 + t111 * t120;
t96 = qJD(2) * qJ(4) + t98;
t92 = t108 * qJD(3) + t110 * t96;
t119 = qJD(2) * t108;
t118 = qJD(2) * t110;
t117 = qJD(1) * qJD(2);
t97 = t111 * t102 - t109 * t120;
t116 = qJD(4) - t97;
t114 = cos(qJ(5));
t112 = sin(qJ(5));
t107 = t110 * qJD(3);
t100 = (t108 * t114 + t110 * t112) * qJD(2);
t99 = t112 * t119 - t114 * t118;
t95 = -qJD(2) * pkin(3) + t116;
t93 = (-pkin(4) * t110 - pkin(3)) * qJD(2) + t116;
t91 = -t108 * t96 + t107;
t90 = pkin(6) * t118 + t92;
t89 = t107 + (-pkin(6) * qJD(2) - t96) * t108;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t115 * t117, -t113 * t117, t98 ^ 2 / 0.2e1 + t97 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, -t95 * t118, t95 * t119, (-t108 * t91 + t110 * t92) * qJD(2), t92 ^ 2 / 0.2e1 + t91 ^ 2 / 0.2e1 + t95 ^ 2 / 0.2e1, t100 ^ 2 / 0.2e1, -t100 * t99, t100 * qJD(5), -t99 * qJD(5), qJD(5) ^ 2 / 0.2e1, (-t112 * t90 + t114 * t89) * qJD(5) + t93 * t99, -(t112 * t89 + t114 * t90) * qJD(5) + t93 * t100;];
T_reg = t1;
