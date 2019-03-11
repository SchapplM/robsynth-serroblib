% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:01:39
% EndTime: 2019-03-09 10:01:39
% DurationCPUTime: 0.12s
% Computational Cost: add. (373->48), mult. (745->98), div. (0->0), fcn. (423->6), ass. (0->39)
t194 = -pkin(2) - pkin(8);
t182 = qJD(1) ^ 2;
t193 = t182 / 0.2e1;
t181 = cos(qJ(2));
t192 = t181 * t182;
t178 = sin(qJ(4));
t180 = cos(qJ(4));
t190 = qJD(1) * t181;
t168 = t180 * qJD(2) - t178 * t190;
t179 = sin(qJ(2));
t189 = t179 * qJD(1);
t171 = qJD(4) + t189;
t184 = -qJ(3) * t179 - pkin(1);
t162 = (t194 * t181 + t184) * qJD(1);
t188 = pkin(7) * t189 + qJD(3);
t163 = pkin(3) * t189 + t194 * qJD(2) + t188;
t183 = -t178 * t162 + t180 * t163;
t153 = t171 * pkin(4) - t168 * qJ(5) + t183;
t167 = t178 * qJD(2) + t180 * t190;
t191 = t180 * t162 + t178 * t163;
t155 = -t167 * qJ(5) + t191;
t176 = sin(pkin(9));
t177 = cos(pkin(9));
t150 = t176 * t153 + t177 * t155;
t170 = -pkin(7) * t190 - qJD(2) * qJ(3);
t187 = qJD(1) * qJD(2);
t165 = pkin(3) * t190 - t170;
t186 = t179 * t187;
t185 = t181 * t187;
t149 = t177 * t153 - t176 * t155;
t158 = t167 * pkin(4) + qJD(5) + t165;
t169 = -qJD(2) * pkin(2) + t188;
t166 = (-pkin(2) * t181 + t184) * qJD(1);
t157 = -t176 * t167 + t177 * t168;
t156 = t177 * t167 + t176 * t168;
t151 = t156 * pkin(5) - t157 * qJ(6) + t158;
t148 = t171 * qJ(6) + t150;
t147 = -t171 * pkin(5) + qJD(6) - t149;
t1 = [t193, 0, 0, t179 ^ 2 * t193, t179 * t192, t186, t185, qJD(2) ^ 2 / 0.2e1, pkin(1) * t192 - pkin(7) * t186, -t182 * pkin(1) * t179 - pkin(7) * t185 (t169 * t179 - t170 * t181) * qJD(1), t169 * qJD(2) + t166 * t190, -t170 * qJD(2) - t166 * t189, t166 ^ 2 / 0.2e1 + t170 ^ 2 / 0.2e1 + t169 ^ 2 / 0.2e1, t168 ^ 2 / 0.2e1, -t168 * t167, t168 * t171, -t167 * t171, t171 ^ 2 / 0.2e1, t165 * t167 + t171 * t183, t165 * t168 - t191 * t171, -t149 * t157 - t150 * t156, t150 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1, -t147 * t171 + t151 * t156, t147 * t157 - t148 * t156, t148 * t171 - t151 * t157, t148 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1 + t147 ^ 2 / 0.2e1;];
T_reg  = t1;
