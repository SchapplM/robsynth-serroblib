% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:07:14
% EndTime: 2019-03-09 05:07:14
% DurationCPUTime: 0.14s
% Computational Cost: add. (236->45), mult. (500->97), div. (0->0), fcn. (299->8), ass. (0->39)
t187 = pkin(4) + pkin(5);
t173 = qJD(1) ^ 2;
t186 = t173 / 0.2e1;
t165 = sin(pkin(10));
t156 = (pkin(1) * t165 + pkin(7)) * qJD(1);
t169 = sin(qJ(3));
t172 = cos(qJ(3));
t183 = t169 * qJD(2) + t172 * t156;
t149 = qJD(3) * pkin(8) + t183;
t166 = cos(pkin(10));
t179 = -pkin(1) * t166 - pkin(2);
t150 = (-pkin(3) * t172 - pkin(8) * t169 + t179) * qJD(1);
t168 = sin(qJ(4));
t171 = cos(qJ(4));
t185 = t171 * t149 + t168 * t150;
t184 = t172 * qJD(2) - t169 * t156;
t182 = qJD(1) * t169;
t181 = t172 * qJD(1);
t180 = qJD(1) * qJD(3);
t160 = -qJD(4) + t181;
t141 = -t160 * qJ(5) + t185;
t178 = -t168 * t149 + t171 * t150;
t177 = qJD(3) * pkin(3) + t184;
t176 = qJD(5) - t178;
t155 = t168 * qJD(3) + t171 * t182;
t175 = t155 * qJ(5) + t177;
t170 = cos(qJ(6));
t167 = sin(qJ(6));
t159 = qJD(6) + t160;
t157 = t179 * qJD(1);
t154 = -t171 * qJD(3) + t168 * t182;
t144 = t167 * t154 + t170 * t155;
t143 = -t170 * t154 + t167 * t155;
t142 = t154 * pkin(4) - t175;
t140 = t160 * pkin(4) + t176;
t139 = -t187 * t154 + t175;
t138 = t154 * pkin(9) + t141;
t137 = -t155 * pkin(9) + t187 * t160 + t176;
t1 = [t186, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t165 ^ 2 / 0.2e1 + t166 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t173, t169 ^ 2 * t186, t169 * t173 * t172, t169 * t180, t172 * t180, qJD(3) ^ 2 / 0.2e1, t184 * qJD(3) - t157 * t181, -t183 * qJD(3) + t157 * t182, t155 ^ 2 / 0.2e1, -t155 * t154, -t155 * t160, t154 * t160, t160 ^ 2 / 0.2e1, -t154 * t177 - t178 * t160, -t155 * t177 + t185 * t160, t140 * t160 + t142 * t154, t140 * t155 - t141 * t154, -t141 * t160 - t142 * t155, t141 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1, t144 ^ 2 / 0.2e1, -t144 * t143, t144 * t159, -t143 * t159, t159 ^ 2 / 0.2e1 (t170 * t137 - t167 * t138) * t159 + t139 * t143 -(t167 * t137 + t170 * t138) * t159 + t139 * t144;];
T_reg  = t1;
