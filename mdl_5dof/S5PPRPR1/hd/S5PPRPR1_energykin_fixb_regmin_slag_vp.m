% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:19
% EndTime: 2019-12-05 15:01:19
% DurationCPUTime: 0.06s
% Computational Cost: add. (65->27), mult. (174->63), div. (0->0), fcn. (120->8), ass. (0->25)
t103 = sin(pkin(8));
t105 = cos(pkin(8));
t115 = qJD(1) * cos(qJ(3));
t116 = qJD(1) * sin(qJ(3));
t117 = t103 * t115 + t105 * t116;
t102 = sin(pkin(9));
t104 = cos(pkin(9));
t93 = qJD(3) * qJ(4) + t117;
t89 = t102 * qJD(2) + t104 * t93;
t114 = qJD(3) * t102;
t113 = qJD(3) * t104;
t112 = -t103 * t116 + t105 * t115;
t111 = qJD(4) - t112;
t110 = qJD(1) ^ 2;
t108 = cos(qJ(5));
t106 = sin(qJ(5));
t101 = t104 * qJD(2);
t95 = (t102 * t108 + t104 * t106) * qJD(3);
t94 = t106 * t114 - t108 * t113;
t92 = -qJD(3) * pkin(3) + t111;
t90 = (-pkin(4) * t104 - pkin(3)) * qJD(3) + t111;
t88 = -t102 * t93 + t101;
t87 = pkin(6) * t113 + t89;
t86 = t101 + (-pkin(6) * qJD(3) - t93) * t102;
t1 = [t110 / 0.2e1, qJD(2) ^ 2 / 0.2e1 + (t103 ^ 2 / 0.2e1 + t105 ^ 2 / 0.2e1) * t110, qJD(3) ^ 2 / 0.2e1, t112 * qJD(3), -t117 * qJD(3), -t92 * t113, t92 * t114, (-t102 * t88 + t104 * t89) * qJD(3), t89 ^ 2 / 0.2e1 + t88 ^ 2 / 0.2e1 + t92 ^ 2 / 0.2e1, t95 ^ 2 / 0.2e1, -t95 * t94, t95 * qJD(5), -t94 * qJD(5), qJD(5) ^ 2 / 0.2e1, (-t106 * t87 + t108 * t86) * qJD(5) + t90 * t94, -(t106 * t86 + t108 * t87) * qJD(5) + t90 * t95;];
T_reg = t1;
