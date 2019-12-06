% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:15
% EndTime: 2019-12-05 15:11:15
% DurationCPUTime: 0.05s
% Computational Cost: add. (61->22), mult. (151->58), div. (0->0), fcn. (88->6), ass. (0->23)
t111 = qJD(3) ^ 2;
t122 = t111 / 0.2e1;
t107 = sin(qJ(4));
t109 = cos(qJ(4));
t108 = sin(qJ(3));
t110 = cos(qJ(3));
t105 = sin(pkin(8));
t119 = qJD(1) * t105;
t115 = t110 * qJD(2) - t108 * t119;
t98 = (-pkin(4) * t109 - qJ(5) * t107 - pkin(3)) * qJD(3) - t115;
t121 = qJD(3) * t98;
t120 = t108 * qJD(2) + t110 * t119;
t106 = cos(pkin(8));
t118 = qJD(1) * t106;
t117 = (-qJD(3) * pkin(3) - t115) * qJD(3);
t116 = qJD(3) * qJD(4);
t101 = qJD(3) * pkin(6) + t120;
t114 = t109 * t101 - t107 * t118;
t113 = -t107 * t101 - t109 * t118;
t112 = qJD(1) ^ 2;
t97 = qJD(4) * qJ(5) + t114;
t96 = -qJD(4) * pkin(4) + qJD(5) - t113;
t1 = [t112 / 0.2e1, qJD(2) ^ 2 / 0.2e1 + (t105 ^ 2 / 0.2e1 + t106 ^ 2 / 0.2e1) * t112, t122, t115 * qJD(3), -t120 * qJD(3), t107 ^ 2 * t122, t107 * t111 * t109, t107 * t116, t109 * t116, qJD(4) ^ 2 / 0.2e1, t113 * qJD(4) - t109 * t117, -t114 * qJD(4) + t107 * t117, -t96 * qJD(4) - t109 * t121, (t107 * t96 + t109 * t97) * qJD(3), t97 * qJD(4) - t107 * t121, t97 ^ 2 / 0.2e1 + t98 ^ 2 / 0.2e1 + t96 ^ 2 / 0.2e1;];
T_reg = t1;
