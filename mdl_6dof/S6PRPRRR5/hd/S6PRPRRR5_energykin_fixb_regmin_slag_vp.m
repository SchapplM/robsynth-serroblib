% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:43:35
% EndTime: 2019-03-08 20:43:35
% DurationCPUTime: 0.08s
% Computational Cost: add. (169->40), mult. (357->88), div. (0->0), fcn. (240->10), ass. (0->39)
t173 = qJD(2) ^ 2;
t187 = t173 / 0.2e1;
t186 = qJD(1) ^ 2 / 0.2e1;
t171 = cos(qJ(4));
t172 = cos(qJ(2));
t183 = qJD(1) * sin(pkin(6));
t176 = -t172 * t183 + qJD(3);
t153 = (-pkin(2) - pkin(8)) * qJD(2) + t176;
t167 = sin(qJ(4));
t164 = cos(pkin(6));
t182 = qJD(1) * t164;
t177 = t171 * t153 - t167 * t182;
t145 = -t171 * qJD(2) * pkin(9) + qJD(4) * pkin(4) + t177;
t181 = qJD(2) * t167;
t184 = t167 * t153 + t171 * t182;
t146 = -pkin(9) * t181 + t184;
t166 = sin(qJ(5));
t170 = cos(qJ(5));
t185 = t166 * t145 + t170 * t146;
t168 = sin(qJ(2));
t157 = qJD(2) * qJ(3) + t168 * t183;
t180 = t157 * qJD(2);
t179 = qJD(2) * qJD(4);
t178 = qJD(2) * t183;
t175 = t170 * t145 - t166 * t146;
t154 = (t166 * t171 + t167 * t170) * qJD(2);
t151 = pkin(4) * t181 + t157;
t169 = cos(qJ(6));
t165 = sin(qJ(6));
t162 = qJD(4) + qJD(5);
t156 = -qJD(2) * pkin(2) + t176;
t155 = (-t166 * t167 + t170 * t171) * qJD(2);
t152 = qJD(6) + t154;
t148 = t169 * t155 + t165 * t162;
t147 = t165 * t155 - t169 * t162;
t142 = t154 * pkin(5) - t155 * pkin(10) + t151;
t141 = t162 * pkin(10) + t185;
t140 = -t162 * pkin(5) - t175;
t1 = [t186, t187, t172 * t178, -t168 * t178, t156 * qJD(2), t180, t164 ^ 2 * t186 + t157 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1, t171 ^ 2 * t187, -t171 * t173 * t167, t171 * t179, -t167 * t179, qJD(4) ^ 2 / 0.2e1, t177 * qJD(4) + t167 * t180, -t184 * qJD(4) + t171 * t180, t155 ^ 2 / 0.2e1, -t155 * t154, t155 * t162, -t154 * t162, t162 ^ 2 / 0.2e1, t151 * t154 + t175 * t162, t151 * t155 - t185 * t162, t148 ^ 2 / 0.2e1, -t148 * t147, t148 * t152, -t147 * t152, t152 ^ 2 / 0.2e1 (-t165 * t141 + t169 * t142) * t152 + t140 * t147 -(t169 * t141 + t165 * t142) * t152 + t140 * t148;];
T_reg  = t1;
