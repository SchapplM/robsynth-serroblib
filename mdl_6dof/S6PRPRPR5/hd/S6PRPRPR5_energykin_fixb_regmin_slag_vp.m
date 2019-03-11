% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRPR5
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
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:45:02
% EndTime: 2019-03-08 19:45:02
% DurationCPUTime: 0.13s
% Computational Cost: add. (219->47), mult. (544->94), div. (0->0), fcn. (397->10), ass. (0->39)
t189 = pkin(4) + pkin(9);
t175 = sin(qJ(2));
t187 = qJD(1) * sin(pkin(6));
t165 = qJD(2) * qJ(3) + t175 * t187;
t171 = cos(pkin(11));
t186 = qJD(1) * cos(pkin(6));
t167 = t171 * t186;
t169 = sin(pkin(11));
t152 = t167 + (-pkin(8) * qJD(2) - t165) * t169;
t157 = t171 * t165 + t169 * t186;
t184 = qJD(2) * t171;
t153 = pkin(8) * t184 + t157;
t174 = sin(qJ(4));
t177 = cos(qJ(4));
t188 = t174 * t152 + t177 * t153;
t185 = qJD(2) * t169;
t183 = qJD(2) * t187;
t182 = t177 * t152 - t174 * t153;
t147 = -qJD(4) * qJ(5) - t188;
t181 = qJD(5) - t182;
t178 = cos(qJ(2));
t180 = -t178 * t187 + qJD(3);
t162 = (t169 * t177 + t171 * t174) * qJD(2);
t159 = (-pkin(3) * t171 - pkin(2)) * qJD(2) + t180;
t179 = -t162 * qJ(5) + t159;
t176 = cos(qJ(6));
t173 = sin(qJ(6));
t164 = -qJD(2) * pkin(2) + t180;
t161 = t174 * t185 - t177 * t184;
t160 = qJD(6) + t162;
t156 = -t169 * t165 + t167;
t155 = t176 * qJD(4) + t173 * t161;
t154 = t173 * qJD(4) - t176 * t161;
t148 = t161 * pkin(4) + t179;
t146 = -qJD(4) * pkin(4) + t181;
t145 = t189 * t161 + t179;
t144 = -t161 * pkin(5) - t147;
t143 = t162 * pkin(5) - t189 * qJD(4) + t181;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t178 * t183, -t175 * t183, -t164 * t184, t164 * t185 (-t156 * t169 + t157 * t171) * qJD(2), t157 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1, t162 ^ 2 / 0.2e1, -t162 * t161, t162 * qJD(4), -t161 * qJD(4), qJD(4) ^ 2 / 0.2e1, t182 * qJD(4) + t159 * t161, -t188 * qJD(4) + t159 * t162, t146 * t162 + t147 * t161, t146 * qJD(4) - t148 * t161, -t147 * qJD(4) - t148 * t162, t148 ^ 2 / 0.2e1 + t147 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1, t155 ^ 2 / 0.2e1, -t155 * t154, t155 * t160, -t154 * t160, t160 ^ 2 / 0.2e1 (t176 * t143 - t173 * t145) * t160 + t144 * t154 -(t173 * t143 + t176 * t145) * t160 + t144 * t155;];
T_reg  = t1;
