% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:13:54
% EndTime: 2019-03-09 10:13:54
% DurationCPUTime: 0.14s
% Computational Cost: add. (420->50), mult. (1036->102), div. (0->0), fcn. (749->8), ass. (0->43)
t194 = qJD(1) * (pkin(7) + qJ(3));
t193 = pkin(4) + pkin(9);
t182 = qJD(1) ^ 2;
t192 = t182 / 0.2e1;
t181 = cos(qJ(2));
t190 = t181 * t182;
t178 = sin(qJ(2));
t170 = qJD(2) * pkin(2) - t178 * t194;
t171 = t181 * t194;
t174 = sin(pkin(10));
t175 = cos(pkin(10));
t161 = t170 * t175 - t171 * t174;
t168 = (t174 * t181 + t175 * t178) * qJD(1);
t153 = qJD(2) * pkin(3) - pkin(8) * t168 + t161;
t162 = t170 * t174 + t171 * t175;
t167 = (-t174 * t178 + t175 * t181) * qJD(1);
t154 = pkin(8) * t167 + t162;
t177 = sin(qJ(4));
t180 = cos(qJ(4));
t189 = t153 * t177 + t154 * t180;
t188 = qJD(1) * qJD(2);
t187 = t178 * t188;
t186 = t181 * t188;
t185 = t180 * t153 - t154 * t177;
t173 = qJD(2) + qJD(4);
t148 = -qJ(5) * t173 - t189;
t184 = qJD(5) - t185;
t160 = t167 * t177 + t168 * t180;
t172 = qJD(3) + (-pkin(2) * t181 - pkin(1)) * qJD(1);
t163 = -t167 * pkin(3) + t172;
t183 = -t160 * qJ(5) + t163;
t179 = cos(qJ(6));
t176 = sin(qJ(6));
t159 = -t167 * t180 + t168 * t177;
t158 = qJD(6) + t160;
t156 = t159 * t176 + t173 * t179;
t155 = -t159 * t179 + t173 * t176;
t149 = t159 * pkin(4) + t183;
t147 = -pkin(4) * t173 + t184;
t146 = t159 * t193 + t183;
t145 = -pkin(5) * t159 - t148;
t144 = t160 * pkin(5) - t173 * t193 + t184;
t1 = [t192, 0, 0, t178 ^ 2 * t192, t178 * t190, t187, t186, qJD(2) ^ 2 / 0.2e1, pkin(1) * t190 - pkin(7) * t187, -pkin(1) * t178 * t182 - pkin(7) * t186, -t161 * t168 + t162 * t167, t162 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1 + t172 ^ 2 / 0.2e1, t160 ^ 2 / 0.2e1, -t160 * t159, t160 * t173, -t159 * t173, t173 ^ 2 / 0.2e1, t159 * t163 + t173 * t185, t160 * t163 - t173 * t189, t147 * t160 + t148 * t159, t147 * t173 - t149 * t159, -t148 * t173 - t149 * t160, t149 ^ 2 / 0.2e1 + t148 ^ 2 / 0.2e1 + t147 ^ 2 / 0.2e1, t156 ^ 2 / 0.2e1, -t156 * t155, t156 * t158, -t155 * t158, t158 ^ 2 / 0.2e1 (t144 * t179 - t146 * t176) * t158 + t145 * t155 -(t144 * t176 + t146 * t179) * t158 + t145 * t156;];
T_reg  = t1;
