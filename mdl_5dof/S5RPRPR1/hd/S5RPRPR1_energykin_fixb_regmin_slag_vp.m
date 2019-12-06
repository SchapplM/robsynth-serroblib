% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:42
% EndTime: 2019-12-05 17:47:42
% DurationCPUTime: 0.12s
% Computational Cost: add. (130->32), mult. (286->71), div. (0->0), fcn. (166->6), ass. (0->28)
t114 = qJD(1) ^ 2;
t119 = t114 / 0.2e1;
t108 = sin(pkin(8));
t109 = cos(pkin(8));
t113 = cos(qJ(3));
t103 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t115 = -qJ(4) * qJD(1) + t103;
t97 = qJD(3) * pkin(3) + t115 * t113;
t111 = sin(qJ(3));
t99 = t115 * t111;
t90 = t108 * t97 + t109 * t99;
t118 = t114 * qJ(2);
t117 = qJD(3) * t103;
t116 = qJD(1) * qJD(3);
t102 = qJD(4) + (pkin(3) * t111 + qJ(2)) * qJD(1);
t89 = -t108 * t99 + t109 * t97;
t112 = cos(qJ(5));
t110 = sin(qJ(5));
t106 = qJD(3) + qJD(5);
t105 = -pkin(1) * qJD(1) + qJD(2);
t101 = (-t108 * t111 + t109 * t113) * qJD(1);
t100 = (-t108 * t113 - t109 * t111) * qJD(1);
t93 = -pkin(4) * t100 + t102;
t92 = t100 * t110 + t101 * t112;
t91 = -t112 * t100 + t101 * t110;
t88 = pkin(7) * t100 + t90;
t87 = qJD(3) * pkin(4) - pkin(7) * t101 + t89;
t1 = [t119, 0, 0, t105 * qJD(1), t118, qJ(2) ^ 2 * t119 + t105 ^ 2 / 0.2e1, t113 ^ 2 * t119, -t113 * t114 * t111, t113 * t116, -t111 * t116, qJD(3) ^ 2 / 0.2e1, t111 * t118 + t113 * t117, -t111 * t117 + t113 * t118, t100 * t90 - t101 * t89, t90 ^ 2 / 0.2e1 + t89 ^ 2 / 0.2e1 + t102 ^ 2 / 0.2e1, t92 ^ 2 / 0.2e1, -t92 * t91, t92 * t106, -t91 * t106, t106 ^ 2 / 0.2e1, t93 * t91 + (-t110 * t88 + t112 * t87) * t106, t93 * t92 - (t110 * t87 + t112 * t88) * t106;];
T_reg = t1;
