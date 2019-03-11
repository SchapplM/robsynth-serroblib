% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:59:17
% EndTime: 2019-03-08 23:59:17
% DurationCPUTime: 0.11s
% Computational Cost: add. (244->43), mult. (554->93), div. (0->0), fcn. (403->10), ass. (0->41)
t181 = qJD(2) ^ 2;
t195 = t181 / 0.2e1;
t194 = cos(qJ(5));
t171 = qJD(3) + qJD(4);
t177 = sin(qJ(2));
t190 = qJD(1) * sin(pkin(6));
t165 = qJD(2) * pkin(8) + t177 * t190;
t179 = cos(qJ(3));
t189 = qJD(1) * cos(pkin(6));
t168 = t179 * t189;
t176 = sin(qJ(3));
t156 = qJD(3) * pkin(3) + t168 + (-pkin(9) * qJD(2) - t165) * t176;
t187 = qJD(2) * t179;
t191 = t179 * t165 + t176 * t189;
t157 = pkin(9) * t187 + t191;
t175 = sin(qJ(4));
t178 = cos(qJ(4));
t192 = t175 * t156 + t178 * t157;
t149 = t171 * pkin(10) + t192;
t180 = cos(qJ(2));
t185 = t180 * t190;
t160 = -t185 + (-pkin(3) * t179 - pkin(2)) * qJD(2);
t188 = qJD(2) * t176;
t162 = t175 * t188 - t178 * t187;
t163 = (t175 * t179 + t176 * t178) * qJD(2);
t152 = t162 * pkin(4) - t163 * pkin(10) + t160;
t174 = sin(qJ(5));
t193 = t194 * t149 + t174 * t152;
t186 = qJD(2) * qJD(3);
t184 = qJD(2) * t190;
t183 = -t174 * t149 + t194 * t152;
t182 = t178 * t156 - t175 * t157;
t148 = -t171 * pkin(4) - t182;
t166 = -qJD(2) * pkin(2) - t185;
t161 = qJD(5) + t162;
t159 = t194 * t163 + t174 * t171;
t158 = t174 * t163 - t194 * t171;
t146 = t158 * pkin(5) + qJD(6) + t148;
t145 = -t158 * qJ(6) + t193;
t144 = t161 * pkin(5) - t159 * qJ(6) + t183;
t1 = [qJD(1) ^ 2 / 0.2e1, t195, t180 * t184, -t177 * t184, t176 ^ 2 * t195, t176 * t181 * t179, t176 * t186, t179 * t186, qJD(3) ^ 2 / 0.2e1 (-t176 * t165 + t168) * qJD(3) - t166 * t187, -qJD(3) * t191 + t166 * t188, t163 ^ 2 / 0.2e1, -t163 * t162, t163 * t171, -t162 * t171, t171 ^ 2 / 0.2e1, t160 * t162 + t171 * t182, t160 * t163 - t171 * t192, t159 ^ 2 / 0.2e1, -t159 * t158, t159 * t161, -t158 * t161, t161 ^ 2 / 0.2e1, t148 * t158 + t161 * t183, t148 * t159 - t161 * t193, -t144 * t159 - t145 * t158, t145 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1;];
T_reg  = t1;
