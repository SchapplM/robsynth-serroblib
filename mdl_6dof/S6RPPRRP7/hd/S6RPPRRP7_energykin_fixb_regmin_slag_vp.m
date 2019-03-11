% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:13:53
% EndTime: 2019-03-09 02:13:53
% DurationCPUTime: 0.08s
% Computational Cost: add. (213->39), mult. (458->85), div. (0->0), fcn. (280->6), ass. (0->34)
t142 = qJD(1) ^ 2;
t150 = t142 / 0.2e1;
t149 = cos(qJ(5));
t137 = sin(pkin(9));
t129 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t145 = -pkin(7) * qJD(1) + t129;
t122 = t145 * t137;
t138 = cos(pkin(9));
t123 = t145 * t138;
t140 = sin(qJ(4));
t141 = cos(qJ(4));
t147 = t141 * t122 + t140 * t123;
t115 = qJD(4) * pkin(8) + t147;
t125 = (t137 * t141 + t138 * t140) * qJD(1);
t126 = (-t137 * t140 + t138 * t141) * qJD(1);
t131 = qJD(1) * qJ(2) + qJD(3);
t146 = qJD(1) * t137;
t127 = pkin(3) * t146 + t131;
t116 = t125 * pkin(4) - t126 * pkin(8) + t127;
t139 = sin(qJ(5));
t148 = t149 * t115 + t139 * t116;
t144 = -t139 * t115 + t149 * t116;
t143 = -t140 * t122 + t141 * t123;
t114 = -qJD(4) * pkin(4) - t143;
t135 = t138 ^ 2;
t134 = t137 ^ 2;
t132 = -qJD(1) * pkin(1) + qJD(2);
t124 = qJD(5) + t125;
t118 = t139 * qJD(4) + t149 * t126;
t117 = -t149 * qJD(4) + t139 * t126;
t110 = t117 * pkin(5) + qJD(6) + t114;
t109 = -t117 * qJ(6) + t148;
t108 = t124 * pkin(5) - t118 * qJ(6) + t144;
t1 = [t150, 0, 0, t132 * qJD(1), t142 * qJ(2), qJ(2) ^ 2 * t150 + t132 ^ 2 / 0.2e1, t131 * t146, t131 * t138 * qJD(1) (-t134 - t135) * t129 * qJD(1), t131 ^ 2 / 0.2e1 + (t134 / 0.2e1 + t135 / 0.2e1) * t129 ^ 2, t126 ^ 2 / 0.2e1, -t126 * t125, t126 * qJD(4), -t125 * qJD(4), qJD(4) ^ 2 / 0.2e1, t143 * qJD(4) + t127 * t125, -t147 * qJD(4) + t127 * t126, t118 ^ 2 / 0.2e1, -t118 * t117, t118 * t124, -t117 * t124, t124 ^ 2 / 0.2e1, t114 * t117 + t144 * t124, t114 * t118 - t148 * t124, -t108 * t118 - t109 * t117, t109 ^ 2 / 0.2e1 + t108 ^ 2 / 0.2e1 + t110 ^ 2 / 0.2e1;];
T_reg  = t1;
