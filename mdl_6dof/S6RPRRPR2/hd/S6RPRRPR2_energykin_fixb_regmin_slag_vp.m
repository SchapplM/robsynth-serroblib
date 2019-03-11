% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:02:33
% EndTime: 2019-03-09 05:02:33
% DurationCPUTime: 0.12s
% Computational Cost: add. (330->46), mult. (724->101), div. (0->0), fcn. (484->10), ass. (0->42)
t184 = qJD(1) ^ 2;
t195 = t184 / 0.2e1;
t194 = cos(qJ(4));
t180 = sin(qJ(4));
t181 = sin(qJ(3));
t191 = qJD(1) * t181;
t166 = qJD(3) * t180 + t191 * t194;
t183 = cos(qJ(3));
t190 = t183 * qJD(1);
t171 = -qJD(4) + t190;
t176 = sin(pkin(10));
t167 = (pkin(1) * t176 + pkin(7)) * qJD(1);
t192 = qJD(2) * t181 + t167 * t183;
t161 = qJD(3) * pkin(8) + t192;
t178 = cos(pkin(10));
t188 = -pkin(1) * t178 - pkin(2);
t162 = (-pkin(3) * t183 - pkin(8) * t181 + t188) * qJD(1);
t187 = -t161 * t180 + t162 * t194;
t149 = -pkin(4) * t171 - qJ(5) * t166 + t187;
t165 = -qJD(3) * t194 + t180 * t191;
t193 = t161 * t194 + t162 * t180;
t152 = -qJ(5) * t165 + t193;
t175 = sin(pkin(11));
t177 = cos(pkin(11));
t144 = t149 * t175 + t152 * t177;
t189 = qJD(1) * qJD(3);
t143 = t149 * t177 - t152 * t175;
t186 = qJD(2) * t183 - t167 * t181;
t160 = -qJD(3) * pkin(3) - t186;
t153 = pkin(4) * t165 + qJD(5) + t160;
t182 = cos(qJ(6));
t179 = sin(qJ(6));
t169 = -qJD(6) + t171;
t168 = t188 * qJD(1);
t156 = -t165 * t175 + t166 * t177;
t155 = -t165 * t177 - t166 * t175;
t150 = -pkin(5) * t155 + t153;
t148 = t155 * t179 + t156 * t182;
t147 = -t155 * t182 + t156 * t179;
t142 = pkin(9) * t155 + t144;
t141 = -pkin(5) * t171 - pkin(9) * t156 + t143;
t1 = [t195, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t176 ^ 2 / 0.2e1 + t178 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t184, t181 ^ 2 * t195, t181 * t184 * t183, t181 * t189, t183 * t189, qJD(3) ^ 2 / 0.2e1, qJD(3) * t186 - t168 * t190, -qJD(3) * t192 + t168 * t191, t166 ^ 2 / 0.2e1, -t166 * t165, -t166 * t171, t165 * t171, t171 ^ 2 / 0.2e1, t160 * t165 - t171 * t187, t160 * t166 + t171 * t193, -t143 * t156 + t144 * t155, t144 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1 + t153 ^ 2 / 0.2e1, t148 ^ 2 / 0.2e1, -t148 * t147, -t148 * t169, t147 * t169, t169 ^ 2 / 0.2e1 -(t141 * t182 - t142 * t179) * t169 + t150 * t147 (t141 * t179 + t142 * t182) * t169 + t150 * t148;];
T_reg  = t1;
