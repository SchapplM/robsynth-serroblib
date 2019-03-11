% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:52:28
% EndTime: 2019-03-09 09:52:28
% DurationCPUTime: 0.15s
% Computational Cost: add. (410->47), mult. (947->95), div. (0->0), fcn. (640->6), ass. (0->40)
t191 = -pkin(4) - pkin(5);
t176 = qJD(1) ^ 2;
t190 = t176 / 0.2e1;
t189 = cos(qJ(4));
t188 = pkin(7) + qJ(3);
t175 = cos(qJ(2));
t187 = t175 * t176;
t171 = sin(pkin(9));
t172 = cos(pkin(9));
t184 = qJD(1) * t175;
t174 = sin(qJ(2));
t185 = qJD(1) * t174;
t162 = -t171 * t185 + t172 * t184;
t163 = (t171 * t175 + t172 * t174) * qJD(1);
t168 = qJD(3) + (-pkin(2) * t175 - pkin(1)) * qJD(1);
t149 = -t162 * pkin(3) - t163 * pkin(8) + t168;
t166 = qJD(2) * pkin(2) - t188 * t185;
t167 = t188 * t184;
t155 = t171 * t166 + t172 * t167;
t153 = qJD(2) * pkin(8) + t155;
t173 = sin(qJ(4));
t186 = t173 * t149 + t189 * t153;
t154 = t172 * t166 - t171 * t167;
t183 = qJD(1) * qJD(2);
t161 = qJD(4) - t162;
t146 = t161 * qJ(5) + t186;
t182 = t174 * t183;
t181 = t175 * t183;
t180 = qJD(2) * pkin(3) + t154;
t179 = t189 * t149 - t173 * t153;
t178 = qJD(5) - t179;
t157 = t173 * qJD(2) + t189 * t163;
t177 = t157 * qJ(5) + t180;
t156 = -t189 * qJD(2) + t173 * t163;
t147 = t156 * pkin(4) - t177;
t145 = -t161 * pkin(4) + t178;
t144 = t191 * t156 + qJD(6) + t177;
t143 = t156 * qJ(6) + t146;
t142 = -t157 * qJ(6) + t191 * t161 + t178;
t1 = [t190, 0, 0, t174 ^ 2 * t190, t174 * t187, t182, t181, qJD(2) ^ 2 / 0.2e1, pkin(1) * t187 - pkin(7) * t182, -t176 * pkin(1) * t174 - pkin(7) * t181, -t154 * t163 + t155 * t162, t155 ^ 2 / 0.2e1 + t154 ^ 2 / 0.2e1 + t168 ^ 2 / 0.2e1, t157 ^ 2 / 0.2e1, -t157 * t156, t157 * t161, -t156 * t161, t161 ^ 2 / 0.2e1, -t156 * t180 + t179 * t161, -t157 * t180 - t186 * t161, -t145 * t161 + t147 * t156, t145 * t157 - t146 * t156, t146 * t161 - t147 * t157, t146 ^ 2 / 0.2e1 + t147 ^ 2 / 0.2e1 + t145 ^ 2 / 0.2e1, -t142 * t161 - t144 * t156, t143 * t161 + t144 * t157, -t142 * t157 + t143 * t156, t143 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1;];
T_reg  = t1;
