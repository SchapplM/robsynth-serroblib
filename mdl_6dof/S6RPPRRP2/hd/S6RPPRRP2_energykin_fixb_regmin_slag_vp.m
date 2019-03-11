% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:01:16
% EndTime: 2019-03-09 02:01:16
% DurationCPUTime: 0.09s
% Computational Cost: add. (299->47), mult. (686->96), div. (0->0), fcn. (453->8), ass. (0->36)
t156 = sin(pkin(9));
t150 = (pkin(1) * t156 + qJ(3)) * qJD(1);
t157 = cos(pkin(10));
t153 = t157 * qJD(2);
t155 = sin(pkin(10));
t138 = t153 + (-pkin(7) * qJD(1) - t150) * t155;
t143 = t155 * qJD(2) + t157 * t150;
t168 = qJD(1) * t157;
t139 = pkin(7) * t168 + t143;
t160 = sin(qJ(4));
t162 = cos(qJ(4));
t170 = t160 * t138 + t162 * t139;
t132 = qJD(4) * pkin(8) + t170;
t158 = cos(pkin(9));
t167 = -pkin(1) * t158 - pkin(2);
t145 = qJD(3) + (-pkin(3) * t157 + t167) * qJD(1);
t169 = qJD(1) * t155;
t146 = t160 * t169 - t162 * t168;
t147 = (t155 * t162 + t157 * t160) * qJD(1);
t134 = t146 * pkin(4) - t147 * pkin(8) + t145;
t159 = sin(qJ(5));
t161 = cos(qJ(5));
t171 = t161 * t132 + t159 * t134;
t166 = t162 * t138 - t160 * t139;
t165 = -t159 * t132 + t161 * t134;
t131 = -qJD(4) * pkin(4) - t166;
t163 = qJD(1) ^ 2;
t149 = t167 * qJD(1) + qJD(3);
t144 = qJD(5) + t146;
t142 = -t155 * t150 + t153;
t141 = t159 * qJD(4) + t161 * t147;
t140 = -t161 * qJD(4) + t159 * t147;
t129 = t140 * pkin(5) - t141 * qJ(6) + t131;
t128 = t144 * qJ(6) + t171;
t127 = -t144 * pkin(5) + qJD(6) - t165;
t1 = [t163 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t156 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t163, -t149 * t168, t149 * t169 (-t142 * t155 + t143 * t157) * qJD(1), t143 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1, t147 ^ 2 / 0.2e1, -t147 * t146, t147 * qJD(4), -t146 * qJD(4), qJD(4) ^ 2 / 0.2e1, t166 * qJD(4) + t145 * t146, -t170 * qJD(4) + t145 * t147, t141 ^ 2 / 0.2e1, -t141 * t140, t141 * t144, -t140 * t144, t144 ^ 2 / 0.2e1, t131 * t140 + t165 * t144, t131 * t141 - t171 * t144, -t127 * t144 + t129 * t140, t127 * t141 - t128 * t140, t128 * t144 - t129 * t141, t128 ^ 2 / 0.2e1 + t129 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1;];
T_reg  = t1;
