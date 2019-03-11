% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:23:43
% EndTime: 2019-03-09 08:23:44
% DurationCPUTime: 0.13s
% Computational Cost: add. (279->54), mult. (628->100), div. (0->0), fcn. (371->6), ass. (0->39)
t173 = sin(pkin(9));
t174 = cos(pkin(9));
t176 = sin(qJ(2));
t187 = qJD(1) * t176;
t161 = -t174 * qJD(2) + t173 * t187;
t193 = (pkin(3) + qJ(5)) * t161;
t192 = pkin(4) + pkin(8);
t179 = qJD(1) ^ 2;
t191 = t179 / 0.2e1;
t189 = -pkin(5) - qJ(4);
t178 = cos(qJ(2));
t188 = t178 * t179;
t160 = (-pkin(2) * t178 - qJ(3) * t176 - pkin(1)) * qJD(1);
t186 = t178 * qJD(1);
t166 = pkin(7) * t186 + qJD(2) * qJ(3);
t157 = t173 * t160 + t174 * t166;
t185 = qJD(1) * qJD(2);
t184 = t176 * t185;
t183 = t178 * t185;
t156 = t174 * t160 - t173 * t166;
t182 = qJD(2) * pkin(2) - pkin(7) * t187 - qJD(3);
t153 = qJ(4) * t186 - t157;
t152 = pkin(3) * t186 + qJD(4) - t156;
t181 = qJ(5) * t186 + t152;
t162 = t173 * qJD(2) + t174 * t187;
t180 = t162 * qJ(4) + t182;
t177 = cos(qJ(6));
t175 = sin(qJ(6));
t167 = -qJD(6) + t186;
t155 = t177 * t161 + t175 * t162;
t154 = t175 * t161 - t177 * t162;
t151 = t161 * pkin(3) - t180;
t150 = -t161 * pkin(4) + qJD(5) - t153;
t149 = t180 - t193;
t148 = t162 * pkin(4) + t181;
t147 = t189 * t162 - t182 + t193;
t146 = -t192 * t161 + t189 * t186 + qJD(5) + t157;
t145 = t192 * t162 + t181;
t1 = [t191, 0, 0, t176 ^ 2 * t191, t176 * t188, t184, t183, qJD(2) ^ 2 / 0.2e1, pkin(1) * t188 - pkin(7) * t184, -t179 * pkin(1) * t176 - pkin(7) * t183, -t156 * t186 - t161 * t182, t157 * t186 - t162 * t182, -t156 * t162 - t157 * t161, t157 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1 + t182 ^ 2 / 0.2e1, t152 * t162 + t153 * t161, -t151 * t161 - t152 * t186, -t151 * t162 + t153 * t186, t151 ^ 2 / 0.2e1 + t153 ^ 2 / 0.2e1 + t152 ^ 2 / 0.2e1, t149 * t162 - t150 * t186, -t148 * t162 + t150 * t161, t148 * t186 - t149 * t161, t148 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1 + t150 ^ 2 / 0.2e1, t155 ^ 2 / 0.2e1, -t155 * t154, -t155 * t167, t154 * t167, t167 ^ 2 / 0.2e1 -(-t175 * t145 + t177 * t146) * t167 + t147 * t154 (t177 * t145 + t175 * t146) * t167 + t147 * t155;];
T_reg  = t1;
