% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:27
% EndTime: 2019-12-31 20:07:27
% DurationCPUTime: 0.12s
% Computational Cost: add. (259->39), mult. (597->84), div. (0->0), fcn. (377->6), ass. (0->33)
t146 = qJD(1) ^ 2;
t155 = t146 / 0.2e1;
t145 = cos(qJ(2));
t154 = t145 * t146;
t143 = sin(qJ(2));
t130 = (-pkin(2) * t145 - qJ(3) * t143 - pkin(1)) * qJD(1);
t151 = t145 * qJD(1);
t135 = pkin(6) * t151 + qJD(2) * qJ(3);
t140 = sin(pkin(8));
t141 = cos(pkin(8));
t124 = t141 * t130 - t140 * t135;
t152 = qJD(1) * t143;
t132 = t140 * qJD(2) + t141 * t152;
t119 = -pkin(3) * t151 - t132 * pkin(7) + t124;
t125 = t140 * t130 + t141 * t135;
t131 = -t141 * qJD(2) + t140 * t152;
t121 = -t131 * pkin(7) + t125;
t142 = sin(qJ(4));
t144 = cos(qJ(4));
t153 = t142 * t119 + t144 * t121;
t150 = qJD(1) * qJD(2);
t149 = t143 * t150;
t148 = t145 * t150;
t134 = -qJD(2) * pkin(2) + pkin(6) * t152 + qJD(3);
t147 = t144 * t119 - t142 * t121;
t126 = t131 * pkin(3) + t134;
t136 = -qJD(4) + t151;
t123 = -t142 * t131 + t144 * t132;
t122 = t144 * t131 + t142 * t132;
t117 = t122 * pkin(4) - t123 * qJ(5) + t126;
t116 = -t136 * qJ(5) + t153;
t115 = t136 * pkin(4) + qJD(5) - t147;
t1 = [t155, 0, 0, t143 ^ 2 * t155, t143 * t154, t149, t148, qJD(2) ^ 2 / 0.2e1, pkin(1) * t154 - pkin(6) * t149, -t146 * pkin(1) * t143 - pkin(6) * t148, -t124 * t151 + t134 * t131, t125 * t151 + t134 * t132, -t124 * t132 - t125 * t131, t125 ^ 2 / 0.2e1 + t124 ^ 2 / 0.2e1 + t134 ^ 2 / 0.2e1, t123 ^ 2 / 0.2e1, -t123 * t122, -t123 * t136, t122 * t136, t136 ^ 2 / 0.2e1, t126 * t122 - t147 * t136, t126 * t123 + t153 * t136, t115 * t136 + t117 * t122, t115 * t123 - t116 * t122, -t116 * t136 - t117 * t123, t116 ^ 2 / 0.2e1 + t117 ^ 2 / 0.2e1 + t115 ^ 2 / 0.2e1;];
T_reg = t1;
