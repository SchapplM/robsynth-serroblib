% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:09:36
% EndTime: 2019-03-09 03:09:36
% DurationCPUTime: 0.14s
% Computational Cost: add. (349->45), mult. (755->98), div. (0->0), fcn. (467->8), ass. (0->37)
t169 = qJD(1) ^ 2;
t179 = t169 / 0.2e1;
t162 = sin(pkin(9));
t155 = (pkin(1) * t162 + pkin(7)) * qJD(1);
t166 = sin(qJ(3));
t168 = cos(qJ(3));
t177 = t166 * qJD(2) + t168 * t155;
t148 = qJD(3) * qJ(4) + t177;
t164 = cos(pkin(9));
t173 = -pkin(1) * t164 - pkin(2);
t149 = (-pkin(3) * t168 - qJ(4) * t166 + t173) * qJD(1);
t161 = sin(pkin(10));
t163 = cos(pkin(10));
t139 = -t148 * t161 + t163 * t149;
t176 = qJD(1) * t166;
t152 = qJD(3) * t161 + t163 * t176;
t175 = qJD(1) * t168;
t136 = -pkin(4) * t175 - pkin(8) * t152 + t139;
t140 = t163 * t148 + t161 * t149;
t151 = -t163 * qJD(3) + t161 * t176;
t138 = -pkin(8) * t151 + t140;
t165 = sin(qJ(5));
t167 = cos(qJ(5));
t178 = t165 * t136 + t167 * t138;
t174 = qJD(1) * qJD(3);
t172 = qJD(2) * t168 - t166 * t155;
t171 = t136 * t167 - t165 * t138;
t147 = -qJD(3) * pkin(3) + qJD(4) - t172;
t141 = pkin(4) * t151 + t147;
t157 = -qJD(5) + t175;
t156 = t173 * qJD(1);
t143 = -t151 * t165 + t152 * t167;
t142 = t167 * t151 + t152 * t165;
t134 = pkin(5) * t142 - qJ(6) * t143 + t141;
t133 = -qJ(6) * t157 + t178;
t132 = pkin(5) * t157 + qJD(6) - t171;
t1 = [t179, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t162 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t169, t166 ^ 2 * t179, t166 * t169 * t168, t166 * t174, t168 * t174, qJD(3) ^ 2 / 0.2e1, t172 * qJD(3) - t156 * t175, -t177 * qJD(3) + t156 * t176, -t139 * t175 + t147 * t151, t140 * t175 + t147 * t152, -t139 * t152 - t140 * t151, t140 ^ 2 / 0.2e1 + t139 ^ 2 / 0.2e1 + t147 ^ 2 / 0.2e1, t143 ^ 2 / 0.2e1, -t143 * t142, -t143 * t157, t142 * t157, t157 ^ 2 / 0.2e1, t141 * t142 - t157 * t171, t141 * t143 + t178 * t157, t132 * t157 + t134 * t142, t132 * t143 - t133 * t142, -t133 * t157 - t134 * t143, t133 ^ 2 / 0.2e1 + t134 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1;];
T_reg  = t1;
