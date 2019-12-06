% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRP5
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
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:40
% EndTime: 2019-12-05 15:38:40
% DurationCPUTime: 0.06s
% Computational Cost: add. (109->29), mult. (261->66), div. (0->0), fcn. (157->6), ass. (0->26)
t106 = sin(qJ(4));
t108 = cos(qJ(4));
t104 = sin(pkin(8));
t107 = sin(qJ(2));
t100 = qJD(2) * qJ(3) + t107 * qJD(1);
t112 = pkin(6) * qJD(2) + t100;
t93 = t112 * t104;
t105 = cos(pkin(8));
t94 = t112 * t105;
t116 = -t106 * t93 + t108 * t94;
t115 = qJD(2) * t104;
t114 = qJD(2) * t105;
t113 = qJD(1) * qJD(2);
t109 = cos(qJ(2));
t111 = -t109 * qJD(1) + qJD(3);
t110 = -t106 * t94 - t108 * t93;
t97 = (-pkin(3) * t105 - pkin(2)) * qJD(2) + t111;
t103 = t105 ^ 2;
t102 = t104 ^ 2;
t98 = -qJD(2) * pkin(2) + t111;
t96 = (t104 * t108 + t105 * t106) * qJD(2);
t95 = t106 * t115 - t108 * t114;
t90 = qJD(4) * qJ(5) + t116;
t89 = -qJD(4) * pkin(4) + qJD(5) - t110;
t88 = t95 * pkin(4) - t96 * qJ(5) + t97;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t109 * t113, -t107 * t113, -t98 * t114, t98 * t115, (t102 + t103) * t100 * qJD(2), t98 ^ 2 / 0.2e1 + (t103 / 0.2e1 + t102 / 0.2e1) * t100 ^ 2, t96 ^ 2 / 0.2e1, -t96 * t95, t96 * qJD(4), -t95 * qJD(4), qJD(4) ^ 2 / 0.2e1, t110 * qJD(4) + t97 * t95, -t116 * qJD(4) + t97 * t96, -t89 * qJD(4) + t88 * t95, t89 * t96 - t90 * t95, t90 * qJD(4) - t88 * t96, t90 ^ 2 / 0.2e1 + t88 ^ 2 / 0.2e1 + t89 ^ 2 / 0.2e1;];
T_reg = t1;
