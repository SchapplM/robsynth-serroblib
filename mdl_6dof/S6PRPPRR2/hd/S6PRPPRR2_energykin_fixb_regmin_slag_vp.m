% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPPRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:20:00
% EndTime: 2019-03-08 19:20:00
% DurationCPUTime: 0.06s
% Computational Cost: add. (108->34), mult. (248->73), div. (0->0), fcn. (161->10), ass. (0->35)
t158 = qJD(2) ^ 2;
t169 = t158 / 0.2e1;
t157 = cos(qJ(2));
t167 = qJD(1) * sin(pkin(6));
t142 = qJD(2) * pkin(2) + t157 * t167;
t149 = sin(pkin(11));
t151 = cos(pkin(11));
t154 = sin(qJ(2));
t162 = t154 * t167;
t138 = t151 * t142 - t149 * t162;
t160 = qJD(4) - t138;
t134 = (-pkin(3) - pkin(8)) * qJD(2) + t160;
t146 = cos(pkin(6)) * qJD(1) + qJD(3);
t153 = sin(qJ(5));
t156 = cos(qJ(5));
t168 = t153 * t134 + t156 * t146;
t166 = qJD(2) * t156;
t139 = t149 * t142 + t151 * t162;
t136 = qJD(2) * qJ(4) + t139;
t165 = t136 * qJD(2);
t164 = t153 * qJD(2);
t163 = qJD(2) * qJD(5);
t161 = qJD(2) * t167;
t159 = t156 * t134 - t153 * t146;
t155 = cos(qJ(6));
t152 = sin(qJ(6));
t147 = qJD(6) + t164;
t145 = t146 ^ 2 / 0.2e1;
t141 = t152 * qJD(5) + t155 * t166;
t140 = -t155 * qJD(5) + t152 * t166;
t135 = -qJD(2) * pkin(3) + t160;
t132 = (pkin(5) * t153 - pkin(9) * t156 + qJ(4)) * qJD(2) + t139;
t131 = qJD(5) * pkin(9) + t168;
t130 = -qJD(5) * pkin(5) - t159;
t1 = [qJD(1) ^ 2 / 0.2e1, t169, t157 * t161, -t154 * t161, t139 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1 + t145, t135 * qJD(2), t165, t145 + t136 ^ 2 / 0.2e1 + t135 ^ 2 / 0.2e1, t156 ^ 2 * t169, -t156 * t158 * t153, t156 * t163, -t153 * t163, qJD(5) ^ 2 / 0.2e1, t159 * qJD(5) + t136 * t164, -t168 * qJD(5) + t156 * t165, t141 ^ 2 / 0.2e1, -t141 * t140, t141 * t147, -t140 * t147, t147 ^ 2 / 0.2e1 (-t152 * t131 + t155 * t132) * t147 + t130 * t140 -(t155 * t131 + t152 * t132) * t147 + t130 * t141;];
T_reg  = t1;
