% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:16
% EndTime: 2019-12-05 16:52:16
% DurationCPUTime: 0.07s
% Computational Cost: add. (117->28), mult. (264->69), div. (0->0), fcn. (157->6), ass. (0->29)
t114 = qJD(2) ^ 2;
t124 = t114 / 0.2e1;
t112 = cos(qJ(3));
t110 = sin(qJ(2));
t104 = qJD(2) * pkin(6) + t110 * qJD(1);
t116 = pkin(7) * qJD(2) + t104;
t100 = t116 * t112;
t108 = sin(qJ(4));
t111 = cos(qJ(4));
t109 = sin(qJ(3));
t99 = qJD(3) * pkin(3) - t116 * t109;
t123 = t111 * t100 + t108 * t99;
t122 = qJD(2) * t109;
t121 = qJD(2) * t112;
t120 = qJD(3) * t104;
t113 = cos(qJ(2));
t119 = t113 * qJD(1);
t118 = qJD(1) * qJD(2);
t117 = qJD(2) * qJD(3);
t115 = -t108 * t100 + t111 * t99;
t103 = -t119 + (-pkin(3) * t112 - pkin(2)) * qJD(2);
t107 = qJD(3) + qJD(4);
t105 = -qJD(2) * pkin(2) - t119;
t102 = (t108 * t112 + t109 * t111) * qJD(2);
t101 = t108 * t122 - t111 * t121;
t96 = t101 * pkin(4) - t102 * qJ(5) + t103;
t95 = t107 * qJ(5) + t123;
t94 = -t107 * pkin(4) + qJD(5) - t115;
t1 = [qJD(1) ^ 2 / 0.2e1, t124, t113 * t118, -t110 * t118, t109 ^ 2 * t124, t109 * t114 * t112, t109 * t117, t112 * t117, qJD(3) ^ 2 / 0.2e1, -t105 * t121 - t109 * t120, t105 * t122 - t112 * t120, t102 ^ 2 / 0.2e1, -t102 * t101, t102 * t107, -t101 * t107, t107 ^ 2 / 0.2e1, t103 * t101 + t115 * t107, t103 * t102 - t123 * t107, t96 * t101 - t94 * t107, -t95 * t101 + t94 * t102, -t96 * t102 + t95 * t107, t95 ^ 2 / 0.2e1 + t96 ^ 2 / 0.2e1 + t94 ^ 2 / 0.2e1;];
T_reg = t1;
