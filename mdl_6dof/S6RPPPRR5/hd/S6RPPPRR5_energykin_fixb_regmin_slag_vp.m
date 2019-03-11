% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPPRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:41
% EndTime: 2019-03-09 01:37:41
% DurationCPUTime: 0.06s
% Computational Cost: add. (126->33), mult. (211->75), div. (0->0), fcn. (81->6), ass. (0->28)
t113 = qJD(1) ^ 2;
t119 = t113 / 0.2e1;
t107 = sin(pkin(9));
t108 = cos(pkin(9));
t102 = qJD(1) * qJ(2) + qJD(3);
t98 = qJD(1) * pkin(3) + t102;
t99 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t92 = t107 * t98 + t108 * t99;
t110 = sin(qJ(5));
t112 = cos(qJ(5));
t90 = qJD(1) * pkin(7) + t92;
t118 = t110 * qJD(4) + t112 * t90;
t117 = qJD(1) * t110;
t116 = t112 * qJD(1);
t115 = qJD(1) * qJD(5);
t91 = -t107 * t99 + t108 * t98;
t114 = t112 * qJD(4) - t110 * t90;
t111 = cos(qJ(6));
t109 = sin(qJ(6));
t103 = -qJD(1) * pkin(1) + qJD(2);
t100 = -qJD(6) + t116;
t94 = t109 * qJD(5) + t111 * t117;
t93 = -t111 * qJD(5) + t109 * t117;
t89 = -qJD(1) * pkin(4) - t91;
t87 = (-pkin(5) * t112 - pkin(8) * t110 - pkin(4)) * qJD(1) - t91;
t86 = qJD(5) * pkin(8) + t118;
t85 = -qJD(5) * pkin(5) - t114;
t1 = [t119, 0, 0, t103 * qJD(1), t113 * qJ(2), qJ(2) ^ 2 * t119 + t103 ^ 2 / 0.2e1, t102 * qJD(1), -t99 * qJD(1), t99 ^ 2 / 0.2e1 + t102 ^ 2 / 0.2e1, t91 * qJD(1), -t92 * qJD(1), t92 ^ 2 / 0.2e1 + t91 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1, t110 ^ 2 * t119, t110 * t113 * t112, t110 * t115, t112 * t115, qJD(5) ^ 2 / 0.2e1, t114 * qJD(5) - t89 * t116, -t118 * qJD(5) + t89 * t117, t94 ^ 2 / 0.2e1, -t94 * t93, -t94 * t100, t93 * t100, t100 ^ 2 / 0.2e1 -(-t109 * t86 + t111 * t87) * t100 + t85 * t93 (t109 * t87 + t111 * t86) * t100 + t85 * t94;];
T_reg  = t1;
