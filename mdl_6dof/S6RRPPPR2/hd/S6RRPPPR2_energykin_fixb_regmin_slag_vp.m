% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:12:06
% EndTime: 2019-03-09 08:12:06
% DurationCPUTime: 0.13s
% Computational Cost: add. (391->51), mult. (923->103), div. (0->0), fcn. (619->8), ass. (0->44)
t191 = qJD(1) ^ 2;
t203 = t191 / 0.2e1;
t202 = pkin(3) + qJ(5);
t201 = pkin(7) + qJ(3);
t200 = cos(pkin(10));
t190 = cos(qJ(2));
t199 = t190 * t191;
t185 = sin(pkin(9));
t186 = cos(pkin(9));
t197 = qJD(1) * t190;
t188 = sin(qJ(2));
t198 = qJD(1) * t188;
t175 = t185 * t198 - t186 * t197;
t176 = (t185 * t190 + t186 * t188) * qJD(1);
t181 = qJD(3) + (-pkin(2) * t190 - pkin(1)) * qJD(1);
t192 = -t176 * qJ(4) + t181;
t159 = t202 * t175 + t192;
t179 = qJD(2) * pkin(2) - t201 * t198;
t180 = t201 * t197;
t167 = t186 * t179 - t185 * t180;
t193 = qJD(4) - t167;
t160 = t176 * pkin(4) - t202 * qJD(2) + t193;
t184 = sin(pkin(10));
t154 = t200 * t159 + t184 * t160;
t168 = t185 * t179 + t186 * t180;
t196 = qJD(1) * qJD(2);
t166 = -qJD(2) * qJ(4) - t168;
t195 = t188 * t196;
t194 = t190 * t196;
t153 = -t184 * t159 + t200 * t160;
t163 = -t175 * pkin(4) + qJD(5) - t166;
t189 = cos(qJ(6));
t187 = sin(qJ(6));
t174 = qJD(6) + t176;
t171 = t200 * qJD(2) + t184 * t175;
t170 = t184 * qJD(2) - t200 * t175;
t165 = -qJD(2) * pkin(3) + t193;
t164 = t175 * pkin(3) + t192;
t162 = -t187 * t170 + t189 * t171;
t161 = t189 * t170 + t187 * t171;
t155 = t170 * pkin(5) + t163;
t152 = -t170 * pkin(8) + t154;
t151 = t176 * pkin(5) - t171 * pkin(8) + t153;
t1 = [t203, 0, 0, t188 ^ 2 * t203, t188 * t199, t195, t194, qJD(2) ^ 2 / 0.2e1, pkin(1) * t199 - pkin(7) * t195, -t191 * pkin(1) * t188 - pkin(7) * t194, -t167 * t176 - t168 * t175, t168 ^ 2 / 0.2e1 + t167 ^ 2 / 0.2e1 + t181 ^ 2 / 0.2e1, t165 * t176 + t166 * t175, t165 * qJD(2) - t164 * t175, -t166 * qJD(2) - t164 * t176, t164 ^ 2 / 0.2e1 + t166 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1, t153 * t176 + t163 * t170, -t154 * t176 + t163 * t171, -t153 * t171 - t154 * t170, t154 ^ 2 / 0.2e1 + t153 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1, t162 ^ 2 / 0.2e1, -t162 * t161, t162 * t174, -t161 * t174, t174 ^ 2 / 0.2e1 (t189 * t151 - t187 * t152) * t174 + t155 * t161 -(t187 * t151 + t189 * t152) * t174 + t155 * t162;];
T_reg  = t1;
