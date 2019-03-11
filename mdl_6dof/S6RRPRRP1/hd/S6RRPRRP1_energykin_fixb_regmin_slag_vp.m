% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:41:14
% EndTime: 2019-03-09 11:41:14
% DurationCPUTime: 0.12s
% Computational Cost: add. (422->47), mult. (1037->98), div. (0->0), fcn. (764->8), ass. (0->42)
t193 = qJD(1) * (pkin(7) + qJ(3));
t181 = qJD(1) ^ 2;
t192 = t181 / 0.2e1;
t191 = cos(qJ(5));
t180 = cos(qJ(2));
t189 = t180 * t181;
t173 = qJD(2) + qJD(4);
t178 = sin(qJ(2));
t169 = qJD(2) * pkin(2) - t178 * t193;
t170 = t180 * t193;
t174 = sin(pkin(10));
t175 = cos(pkin(10));
t160 = t175 * t169 - t174 * t170;
t167 = (t174 * t180 + t175 * t178) * qJD(1);
t153 = qJD(2) * pkin(3) - t167 * pkin(8) + t160;
t161 = t174 * t169 + t175 * t170;
t166 = (-t174 * t178 + t175 * t180) * qJD(1);
t154 = t166 * pkin(8) + t161;
t177 = sin(qJ(4));
t179 = cos(qJ(4));
t187 = t177 * t153 + t179 * t154;
t146 = t173 * pkin(9) + t187;
t158 = -t179 * t166 + t177 * t167;
t159 = t177 * t166 + t179 * t167;
t171 = qJD(3) + (-pkin(2) * t180 - pkin(1)) * qJD(1);
t162 = -t166 * pkin(3) + t171;
t149 = t158 * pkin(4) - t159 * pkin(9) + t162;
t176 = sin(qJ(5));
t188 = t191 * t146 + t176 * t149;
t186 = qJD(1) * qJD(2);
t185 = t178 * t186;
t184 = t180 * t186;
t183 = -t176 * t146 + t191 * t149;
t182 = t179 * t153 - t177 * t154;
t145 = -t173 * pkin(4) - t182;
t157 = qJD(5) + t158;
t156 = t191 * t159 + t176 * t173;
t155 = t176 * t159 - t191 * t173;
t143 = t155 * pkin(5) + qJD(6) + t145;
t142 = -t155 * qJ(6) + t188;
t141 = t157 * pkin(5) - t156 * qJ(6) + t183;
t1 = [t192, 0, 0, t178 ^ 2 * t192, t178 * t189, t185, t184, qJD(2) ^ 2 / 0.2e1, pkin(1) * t189 - pkin(7) * t185, -t181 * pkin(1) * t178 - pkin(7) * t184, -t160 * t167 + t161 * t166, t161 ^ 2 / 0.2e1 + t160 ^ 2 / 0.2e1 + t171 ^ 2 / 0.2e1, t159 ^ 2 / 0.2e1, -t159 * t158, t159 * t173, -t158 * t173, t173 ^ 2 / 0.2e1, t162 * t158 + t182 * t173, t162 * t159 - t187 * t173, t156 ^ 2 / 0.2e1, -t156 * t155, t156 * t157, -t155 * t157, t157 ^ 2 / 0.2e1, t145 * t155 + t183 * t157, t145 * t156 - t188 * t157, -t141 * t156 - t142 * t155, t142 ^ 2 / 0.2e1 + t141 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1;];
T_reg  = t1;
