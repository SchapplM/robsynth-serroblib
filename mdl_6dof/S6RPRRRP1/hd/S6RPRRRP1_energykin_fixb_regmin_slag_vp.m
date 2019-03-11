% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:57:05
% EndTime: 2019-03-09 05:57:05
% DurationCPUTime: 0.09s
% Computational Cost: add. (329->45), mult. (693->99), div. (0->0), fcn. (449->8), ass. (0->38)
t163 = qJD(1) ^ 2;
t174 = t163 / 0.2e1;
t154 = qJD(3) + qJD(4);
t155 = sin(pkin(10));
t148 = (pkin(1) * t155 + pkin(7)) * qJD(1);
t162 = cos(qJ(3));
t153 = t162 * qJD(2);
t159 = sin(qJ(3));
t139 = qJD(3) * pkin(3) + t153 + (-pkin(8) * qJD(1) - t148) * t159;
t169 = qJD(1) * t162;
t171 = t159 * qJD(2) + t162 * t148;
t142 = pkin(8) * t169 + t171;
t158 = sin(qJ(4));
t161 = cos(qJ(4));
t172 = t158 * t139 + t161 * t142;
t133 = t154 * pkin(9) + t172;
t170 = qJD(1) * t159;
t144 = t158 * t170 - t161 * t169;
t145 = (t158 * t162 + t159 * t161) * qJD(1);
t156 = cos(pkin(10));
t167 = -pkin(1) * t156 - pkin(2);
t146 = (-pkin(3) * t162 + t167) * qJD(1);
t135 = t144 * pkin(4) - t145 * pkin(9) + t146;
t157 = sin(qJ(5));
t160 = cos(qJ(5));
t173 = t160 * t133 + t157 * t135;
t168 = qJD(1) * qJD(3);
t166 = t161 * t139 - t158 * t142;
t165 = -t157 * t133 + t160 * t135;
t132 = -t154 * pkin(4) - t166;
t149 = t167 * qJD(1);
t143 = qJD(5) + t144;
t141 = t160 * t145 + t157 * t154;
t140 = t157 * t145 - t160 * t154;
t130 = t140 * pkin(5) - t141 * qJ(6) + t132;
t129 = t143 * qJ(6) + t173;
t128 = -t143 * pkin(5) + qJD(6) - t165;
t1 = [t174, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t155 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t163, t159 ^ 2 * t174, t159 * t163 * t162, t159 * t168, t162 * t168, qJD(3) ^ 2 / 0.2e1, -t149 * t169 + (-t159 * t148 + t153) * qJD(3), -t171 * qJD(3) + t149 * t170, t145 ^ 2 / 0.2e1, -t145 * t144, t145 * t154, -t144 * t154, t154 ^ 2 / 0.2e1, t146 * t144 + t166 * t154, t146 * t145 - t172 * t154, t141 ^ 2 / 0.2e1, -t141 * t140, t141 * t143, -t140 * t143, t143 ^ 2 / 0.2e1, t132 * t140 + t165 * t143, t132 * t141 - t173 * t143, -t128 * t143 + t130 * t140, t128 * t141 - t129 * t140, t129 * t143 - t130 * t141, t129 ^ 2 / 0.2e1 + t130 ^ 2 / 0.2e1 + t128 ^ 2 / 0.2e1;];
T_reg  = t1;
