% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPPRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:30:13
% EndTime: 2019-03-09 01:30:13
% DurationCPUTime: 0.06s
% Computational Cost: add. (88->34), mult. (182->71), div. (0->0), fcn. (73->6), ass. (0->28)
t111 = qJD(1) ^ 2;
t120 = t111 / 0.2e1;
t108 = sin(qJ(5));
t110 = cos(qJ(5));
t105 = sin(pkin(9));
t98 = (-pkin(1) * t105 - qJ(3)) * qJD(1);
t96 = qJD(4) - t98;
t91 = -qJD(1) * pkin(7) + t96;
t119 = t110 * qJD(2) + t108 * t91;
t118 = qJD(1) * t110;
t117 = t108 * qJD(1);
t116 = qJD(1) * qJD(5);
t106 = cos(pkin(9));
t115 = -pkin(1) * t106 - pkin(2);
t114 = -qJ(4) + t115;
t113 = -t108 * qJD(2) + t110 * t91;
t109 = cos(qJ(6));
t107 = sin(qJ(6));
t103 = qJD(2) ^ 2 / 0.2e1;
t99 = qJD(6) + t117;
t97 = t115 * qJD(1) + qJD(3);
t95 = t107 * qJD(5) + t109 * t118;
t94 = -t109 * qJD(5) + t107 * t118;
t92 = -t114 * qJD(1) - qJD(3);
t89 = -qJD(3) + (pkin(5) * t108 - pkin(8) * t110 - t114) * qJD(1);
t88 = qJD(5) * pkin(8) + t119;
t87 = -qJD(5) * pkin(5) - t113;
t1 = [t120, 0, 0, t103 + (t105 ^ 2 / 0.2e1 + t106 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t111, t97 * qJD(1), -t98 * qJD(1), t103 + t98 ^ 2 / 0.2e1 + t97 ^ 2 / 0.2e1, t96 * qJD(1), t92 * qJD(1), t103 + t92 ^ 2 / 0.2e1 + t96 ^ 2 / 0.2e1, t110 ^ 2 * t120, -t110 * t111 * t108, t110 * t116, -t108 * t116, qJD(5) ^ 2 / 0.2e1, t113 * qJD(5) + t92 * t117, -t119 * qJD(5) + t92 * t118, t95 ^ 2 / 0.2e1, -t95 * t94, t95 * t99, -t94 * t99, t99 ^ 2 / 0.2e1 (-t107 * t88 + t109 * t89) * t99 + t87 * t94 -(t107 * t89 + t109 * t88) * t99 + t87 * t95;];
T_reg  = t1;
