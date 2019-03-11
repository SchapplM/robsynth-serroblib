% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x31]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:09:19
% EndTime: 2019-03-09 04:09:19
% DurationCPUTime: 0.12s
% Computational Cost: add. (319->48), mult. (656->99), div. (0->0), fcn. (432->8), ass. (0->39)
t167 = qJD(1) ^ 2;
t176 = t167 / 0.2e1;
t175 = cos(qJ(5));
t174 = t167 * qJ(2);
t164 = sin(qJ(3));
t166 = cos(qJ(3));
t151 = (pkin(3) * t164 - qJ(4) * t166 + qJ(2)) * qJD(1);
t156 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t152 = qJD(3) * qJ(4) + t164 * t156;
t160 = sin(pkin(10));
t161 = cos(pkin(10));
t140 = t161 * t151 - t160 * t152;
t172 = qJD(1) * t166;
t154 = t160 * qJD(3) + t161 * t172;
t170 = t164 * qJD(1);
t137 = pkin(4) * t170 - t154 * pkin(8) + t140;
t141 = t160 * t151 + t161 * t152;
t153 = -t161 * qJD(3) + t160 * t172;
t139 = -t153 * pkin(8) + t141;
t163 = sin(qJ(5));
t173 = t163 * t137 + t175 * t139;
t171 = qJD(3) * t156;
t169 = qJD(1) * qJD(3);
t168 = t175 * t137 - t163 * t139;
t157 = qJD(5) + t170;
t150 = -qJD(3) * pkin(3) - t166 * t156 + qJD(4);
t145 = t153 * pkin(4) + t150;
t165 = cos(qJ(6));
t162 = sin(qJ(6));
t158 = -qJD(1) * pkin(1) + qJD(2);
t155 = qJD(6) + t157;
t144 = -t163 * t153 + t175 * t154;
t143 = t175 * t153 + t163 * t154;
t134 = t143 * pkin(5) + t145;
t133 = -t162 * t143 + t165 * t144;
t132 = t165 * t143 + t162 * t144;
t131 = -t143 * pkin(9) + t173;
t130 = t157 * pkin(5) - t144 * pkin(9) + t168;
t1 = [t176, 0, 0, t158 * qJD(1), t174, qJ(2) ^ 2 * t176 + t158 ^ 2 / 0.2e1, t166 ^ 2 * t176, -t166 * t167 * t164, t166 * t169, -t164 * t169, qJD(3) ^ 2 / 0.2e1, t164 * t174 + t166 * t171, -t164 * t171 + t166 * t174, t140 * t170 + t150 * t153, -t141 * t170 + t150 * t154, -t140 * t154 - t141 * t153, t141 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1 + t150 ^ 2 / 0.2e1, t144 ^ 2 / 0.2e1, -t144 * t143, t144 * t157, -t143 * t157, t157 ^ 2 / 0.2e1, t145 * t143 + t168 * t157, t145 * t144 - t173 * t157, t133 ^ 2 / 0.2e1, -t133 * t132, t133 * t155, -t132 * t155, t155 ^ 2 / 0.2e1 (t165 * t130 - t162 * t131) * t155 + t134 * t132 -(t162 * t130 + t165 * t131) * t155 + t134 * t133;];
T_reg  = t1;
