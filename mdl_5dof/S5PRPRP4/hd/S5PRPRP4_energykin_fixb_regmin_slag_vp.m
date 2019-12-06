% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:05
% EndTime: 2019-12-05 15:36:05
% DurationCPUTime: 0.05s
% Computational Cost: add. (76->23), mult. (171->59), div. (0->0), fcn. (95->6), ass. (0->24)
t108 = qJD(2) ^ 2;
t116 = t108 / 0.2e1;
t102 = sin(pkin(8));
t103 = cos(pkin(8));
t105 = sin(qJ(2));
t114 = qJD(1) * t105;
t107 = cos(qJ(2));
t98 = qJD(2) * pkin(2) + t107 * qJD(1);
t96 = t102 * t98 + t103 * t114;
t104 = sin(qJ(4));
t106 = cos(qJ(4));
t94 = qJD(2) * pkin(6) + t96;
t115 = t104 * qJD(3) + t106 * t94;
t113 = qJD(2) * t104;
t112 = qJD(2) * t106;
t111 = qJD(1) * qJD(2);
t110 = qJD(2) * qJD(4);
t95 = -t102 * t114 + t103 * t98;
t109 = t106 * qJD(3) - t104 * t94;
t93 = -qJD(2) * pkin(3) - t95;
t91 = (-pkin(4) * t106 - qJ(5) * t104 - pkin(3)) * qJD(2) - t95;
t90 = qJD(4) * qJ(5) + t115;
t89 = -qJD(4) * pkin(4) + qJD(5) - t109;
t1 = [qJD(1) ^ 2 / 0.2e1, t116, t107 * t111, -t105 * t111, t96 ^ 2 / 0.2e1 + t95 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t104 ^ 2 * t116, t104 * t108 * t106, t104 * t110, t106 * t110, qJD(4) ^ 2 / 0.2e1, t109 * qJD(4) - t93 * t112, -t115 * qJD(4) + t93 * t113, -t89 * qJD(4) - t91 * t112, (t104 * t89 + t106 * t90) * qJD(2), t90 * qJD(4) - t91 * t113, t90 ^ 2 / 0.2e1 + t91 ^ 2 / 0.2e1 + t89 ^ 2 / 0.2e1;];
T_reg = t1;
