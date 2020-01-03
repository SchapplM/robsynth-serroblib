% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:12
% EndTime: 2019-12-31 20:18:12
% DurationCPUTime: 0.14s
% Computational Cost: add. (236->39), mult. (610->86), div. (0->0), fcn. (442->8), ass. (0->37)
t162 = qJD(1) * (pkin(6) + qJ(3));
t153 = qJD(1) ^ 2;
t161 = t153 / 0.2e1;
t152 = cos(qJ(2));
t159 = t152 * t153;
t149 = sin(qJ(2));
t140 = qJD(2) * pkin(2) - t149 * t162;
t141 = t152 * t162;
t145 = sin(pkin(9));
t146 = cos(pkin(9));
t131 = t146 * t140 - t145 * t141;
t138 = (t145 * t152 + t146 * t149) * qJD(1);
t124 = qJD(2) * pkin(3) - t138 * pkin(7) + t131;
t132 = t145 * t140 + t146 * t141;
t137 = (-t145 * t149 + t146 * t152) * qJD(1);
t125 = t137 * pkin(7) + t132;
t148 = sin(qJ(4));
t151 = cos(qJ(4));
t158 = t148 * t124 + t151 * t125;
t157 = qJD(1) * qJD(2);
t156 = t149 * t157;
t155 = t152 * t157;
t129 = -t151 * t137 + t148 * t138;
t154 = t151 * t124 - t148 * t125;
t142 = qJD(3) + (-pkin(2) * t152 - pkin(1)) * qJD(1);
t133 = -t137 * pkin(3) + t142;
t150 = cos(qJ(5));
t147 = sin(qJ(5));
t144 = qJD(2) + qJD(4);
t130 = t148 * t137 + t151 * t138;
t128 = qJD(5) + t129;
t127 = t150 * t130 + t147 * t144;
t126 = t147 * t130 - t150 * t144;
t121 = t129 * pkin(4) - t130 * pkin(8) + t133;
t120 = t144 * pkin(8) + t158;
t119 = -t144 * pkin(4) - t154;
t1 = [t161, 0, 0, t149 ^ 2 * t161, t149 * t159, t156, t155, qJD(2) ^ 2 / 0.2e1, pkin(1) * t159 - pkin(6) * t156, -t153 * pkin(1) * t149 - pkin(6) * t155, -t131 * t138 + t132 * t137, t132 ^ 2 / 0.2e1 + t131 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1, t130 ^ 2 / 0.2e1, -t130 * t129, t130 * t144, -t129 * t144, t144 ^ 2 / 0.2e1, t133 * t129 + t154 * t144, t133 * t130 - t158 * t144, t127 ^ 2 / 0.2e1, -t127 * t126, t127 * t128, -t126 * t128, t128 ^ 2 / 0.2e1, (-t147 * t120 + t150 * t121) * t128 + t119 * t126, -(t150 * t120 + t147 * t121) * t128 + t119 * t127;];
T_reg = t1;
