% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:07:44
% EndTime: 2019-03-08 20:07:44
% DurationCPUTime: 0.10s
% Computational Cost: add. (223->44), mult. (553->90), div. (0->0), fcn. (412->10), ass. (0->38)
t189 = cos(qJ(5));
t176 = sin(qJ(2));
t186 = qJD(1) * sin(pkin(6));
t165 = qJD(2) * qJ(3) + t176 * t186;
t172 = cos(pkin(11));
t185 = qJD(1) * cos(pkin(6));
t167 = t172 * t185;
t170 = sin(pkin(11));
t153 = t167 + (-pkin(8) * qJD(2) - t165) * t170;
t158 = t165 * t172 + t170 * t185;
t183 = qJD(2) * t172;
t154 = pkin(8) * t183 + t158;
t175 = sin(qJ(4));
t177 = cos(qJ(4));
t187 = t153 * t175 + t154 * t177;
t146 = qJD(4) * pkin(9) + t187;
t178 = cos(qJ(2));
t179 = -t178 * t186 + qJD(3);
t159 = (-pkin(3) * t172 - pkin(2)) * qJD(2) + t179;
t184 = qJD(2) * t170;
t161 = t175 * t184 - t177 * t183;
t162 = (t170 * t177 + t172 * t175) * qJD(2);
t149 = t161 * pkin(4) - t162 * pkin(9) + t159;
t174 = sin(qJ(5));
t188 = t146 * t189 + t149 * t174;
t182 = qJD(2) * t186;
t181 = -t146 * t174 + t149 * t189;
t180 = t153 * t177 - t154 * t175;
t145 = -qJD(4) * pkin(4) - t180;
t164 = -qJD(2) * pkin(2) + t179;
t160 = qJD(5) + t161;
t157 = -t165 * t170 + t167;
t156 = qJD(4) * t174 + t162 * t189;
t155 = -qJD(4) * t189 + t162 * t174;
t143 = pkin(5) * t155 + qJD(6) + t145;
t142 = -qJ(6) * t155 + t188;
t141 = pkin(5) * t160 - qJ(6) * t156 + t181;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t178 * t182, -t176 * t182, -t164 * t183, t164 * t184 (-t157 * t170 + t158 * t172) * qJD(2), t158 ^ 2 / 0.2e1 + t157 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1, t162 ^ 2 / 0.2e1, -t162 * t161, t162 * qJD(4), -t161 * qJD(4), qJD(4) ^ 2 / 0.2e1, qJD(4) * t180 + t159 * t161, -qJD(4) * t187 + t159 * t162, t156 ^ 2 / 0.2e1, -t156 * t155, t156 * t160, -t155 * t160, t160 ^ 2 / 0.2e1, t145 * t155 + t160 * t181, t145 * t156 - t160 * t188, -t141 * t156 - t142 * t155, t142 ^ 2 / 0.2e1 + t141 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1;];
T_reg  = t1;
