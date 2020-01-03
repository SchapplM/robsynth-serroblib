% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:23
% EndTime: 2020-01-03 11:50:23
% DurationCPUTime: 0.13s
% Computational Cost: add. (146->37), mult. (385->85), div. (0->0), fcn. (236->6), ass. (0->32)
t131 = sin(pkin(8));
t129 = t131 ^ 2;
t137 = qJD(1) ^ 2;
t147 = t129 * t137;
t132 = cos(pkin(8));
t119 = qJD(2) + (-pkin(2) * t132 - pkin(6) * t131 - pkin(1)) * qJD(1);
t136 = cos(qJ(3));
t118 = t136 * t119;
t143 = t132 * qJD(1);
t125 = -qJD(3) + t143;
t134 = sin(qJ(3));
t111 = -t125 * pkin(3) + t118 + (-pkin(7) * t131 * t136 - qJ(2) * t132 * t134) * qJD(1);
t144 = qJD(1) * t131;
t140 = t134 * t144;
t141 = qJ(2) * t143;
t145 = t134 * t119 + t136 * t141;
t114 = -pkin(7) * t140 + t145;
t133 = sin(qJ(4));
t135 = cos(qJ(4));
t146 = t133 * t111 + t135 * t114;
t120 = pkin(3) * t140 + qJ(2) * t144;
t142 = t134 * t147;
t139 = t135 * t111 - t133 * t114;
t130 = t132 ^ 2;
t128 = -qJD(1) * pkin(1) + qJD(2);
t122 = -qJD(4) + t125;
t116 = (-t133 * t134 + t135 * t136) * t144;
t115 = (t133 * t136 + t134 * t135) * t144;
t113 = t115 * pkin(4) + qJD(5) + t120;
t108 = -t115 * qJ(5) + t146;
t107 = -t122 * pkin(4) - t116 * qJ(5) + t139;
t1 = [t137 / 0.2e1, 0, 0, -t128 * t143, t128 * t144, (t129 + t130) * t137 * qJ(2), t128 ^ 2 / 0.2e1 + (t130 / 0.2e1 + t129 / 0.2e1) * qJ(2) ^ 2 * t137, t136 ^ 2 * t147 / 0.2e1, -t136 * t142, -t136 * t125 * t144, t125 * t140, t125 ^ 2 / 0.2e1, qJ(2) * t142 - (-t134 * t141 + t118) * t125, qJ(2) * t136 * t147 + t145 * t125, t116 ^ 2 / 0.2e1, -t116 * t115, -t116 * t122, t115 * t122, t122 ^ 2 / 0.2e1, t120 * t115 - t139 * t122, t120 * t116 + t146 * t122, -t107 * t116 - t108 * t115, t108 ^ 2 / 0.2e1 + t107 ^ 2 / 0.2e1 + t113 ^ 2 / 0.2e1;];
T_reg = t1;
