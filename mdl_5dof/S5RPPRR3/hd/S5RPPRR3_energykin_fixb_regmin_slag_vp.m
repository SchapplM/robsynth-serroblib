% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:28:13
% EndTime: 2020-01-03 11:28:13
% DurationCPUTime: 0.14s
% Computational Cost: add. (127->37), mult. (331->80), div. (0->0), fcn. (215->8), ass. (0->31)
t127 = cos(qJ(4));
t118 = sin(qJ(4));
t114 = sin(pkin(8));
t108 = (pkin(1) * t114 + qJ(3)) * qJD(1);
t115 = cos(pkin(9));
t111 = t115 * qJD(2);
t113 = sin(pkin(9));
t98 = t111 + (-pkin(6) * qJD(1) - t108) * t113;
t101 = t113 * qJD(2) + t115 * t108;
t124 = qJD(1) * t115;
t99 = pkin(6) * t124 + t101;
t126 = t118 * t98 + t127 * t99;
t125 = qJD(1) * t113;
t116 = cos(pkin(8));
t123 = -pkin(1) * t116 - pkin(2);
t122 = -t118 * t99 + t127 * t98;
t103 = qJD(3) + (-pkin(3) * t115 + t123) * qJD(1);
t120 = qJD(1) ^ 2;
t119 = cos(qJ(5));
t117 = sin(qJ(5));
t112 = qJD(4) + qJD(5);
t107 = t123 * qJD(1) + qJD(3);
t105 = (t127 * t113 + t115 * t118) * qJD(1);
t104 = t118 * t125 - t127 * t124;
t100 = -t113 * t108 + t111;
t94 = t104 * pkin(4) + t103;
t93 = -t117 * t104 + t119 * t105;
t92 = t119 * t104 + t117 * t105;
t91 = -t104 * pkin(7) + t126;
t90 = qJD(4) * pkin(4) - t105 * pkin(7) + t122;
t1 = [t120 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t114 ^ 2 / 0.2e1 + t116 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t120, -t107 * t124, t107 * t125, (-t100 * t113 + t101 * t115) * qJD(1), t101 ^ 2 / 0.2e1 + t100 ^ 2 / 0.2e1 + t107 ^ 2 / 0.2e1, t105 ^ 2 / 0.2e1, -t105 * t104, t105 * qJD(4), -t104 * qJD(4), qJD(4) ^ 2 / 0.2e1, t122 * qJD(4) + t103 * t104, -t126 * qJD(4) + t103 * t105, t93 ^ 2 / 0.2e1, -t93 * t92, t93 * t112, -t92 * t112, t112 ^ 2 / 0.2e1, t94 * t92 + (-t117 * t91 + t119 * t90) * t112, t94 * t93 - (t117 * t90 + t119 * t91) * t112;];
T_reg = t1;
