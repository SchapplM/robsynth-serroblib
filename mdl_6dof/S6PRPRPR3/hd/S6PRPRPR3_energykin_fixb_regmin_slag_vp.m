% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:37:38
% EndTime: 2019-03-08 19:37:38
% DurationCPUTime: 0.08s
% Computational Cost: add. (145->39), mult. (339->85), div. (0->0), fcn. (224->10), ass. (0->37)
t186 = -pkin(4) - pkin(9);
t174 = qJD(2) ^ 2;
t185 = t174 / 0.2e1;
t173 = cos(qJ(2));
t183 = qJD(1) * sin(pkin(6));
t159 = qJD(2) * pkin(2) + t173 * t183;
t165 = sin(pkin(11));
t167 = cos(pkin(11));
t170 = sin(qJ(2));
t179 = t170 * t183;
t155 = t165 * t159 + t167 * t179;
t153 = qJD(2) * pkin(8) + t155;
t163 = cos(pkin(6)) * qJD(1) + qJD(3);
t169 = sin(qJ(4));
t172 = cos(qJ(4));
t184 = t172 * t153 + t169 * t163;
t182 = qJD(2) * t172;
t181 = t169 * qJD(2);
t180 = qJD(2) * qJD(4);
t178 = qJD(2) * t183;
t177 = -qJ(5) * t169 - pkin(3);
t176 = -t169 * t153 + t172 * t163;
t154 = t167 * t159 - t165 * t179;
t147 = -qJD(4) * qJ(5) - t184;
t175 = qJD(5) - t176;
t171 = cos(qJ(6));
t168 = sin(qJ(6));
t164 = qJD(6) + t181;
t158 = t171 * qJD(4) - t168 * t182;
t157 = t168 * qJD(4) + t171 * t182;
t152 = -qJD(2) * pkin(3) - t154;
t149 = (-pkin(4) * t172 + t177) * qJD(2) - t154;
t148 = (t186 * t172 + t177) * qJD(2) - t154;
t146 = -qJD(4) * pkin(4) + t175;
t145 = pkin(5) * t182 - t147;
t144 = pkin(5) * t181 + t186 * qJD(4) + t175;
t1 = [qJD(1) ^ 2 / 0.2e1, t185, t173 * t178, -t170 * t178, t155 ^ 2 / 0.2e1 + t154 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1, t169 ^ 2 * t185, t169 * t174 * t172, t169 * t180, t172 * t180, qJD(4) ^ 2 / 0.2e1, t176 * qJD(4) - t152 * t182, -t184 * qJD(4) + t152 * t181 (t146 * t169 - t147 * t172) * qJD(2), t146 * qJD(4) + t149 * t182, -t147 * qJD(4) - t149 * t181, t149 ^ 2 / 0.2e1 + t147 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1, t158 ^ 2 / 0.2e1, -t158 * t157, t158 * t164, -t157 * t164, t164 ^ 2 / 0.2e1 (t171 * t144 - t168 * t148) * t164 + t145 * t157 -(t168 * t144 + t171 * t148) * t164 + t145 * t158;];
T_reg  = t1;
