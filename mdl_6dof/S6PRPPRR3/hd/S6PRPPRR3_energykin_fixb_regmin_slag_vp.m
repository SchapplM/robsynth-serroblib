% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPPRR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:23:02
% EndTime: 2019-03-08 19:23:02
% DurationCPUTime: 0.12s
% Computational Cost: add. (134->34), mult. (267->78), div. (0->0), fcn. (157->10), ass. (0->35)
t158 = qJD(2) ^ 2;
t169 = t158 / 0.2e1;
t168 = qJD(1) ^ 2 / 0.2e1;
t157 = cos(qJ(2));
t166 = qJD(1) * sin(pkin(6));
t161 = -t157 * t166 + qJD(3);
t137 = (-pkin(2) - pkin(3)) * qJD(2) + t161;
t154 = sin(qJ(2));
t143 = qJD(2) * qJ(3) + t154 * t166;
t147 = sin(pkin(11));
t149 = cos(pkin(11));
t135 = t147 * t137 + t149 * t143;
t133 = -qJD(2) * pkin(8) + t135;
t150 = cos(pkin(6));
t145 = -t150 * qJD(1) + qJD(4);
t153 = sin(qJ(5));
t156 = cos(qJ(5));
t167 = t156 * t133 + t153 * t145;
t165 = qJD(2) * t153;
t164 = t156 * qJD(2);
t163 = qJD(2) * qJD(5);
t162 = qJD(2) * t166;
t134 = t149 * t137 - t147 * t143;
t132 = qJD(2) * pkin(4) - t134;
t160 = -t153 * t133 + t156 * t145;
t155 = cos(qJ(6));
t152 = sin(qJ(6));
t146 = qJD(6) + t164;
t142 = t152 * qJD(5) - t155 * t165;
t141 = t155 * qJD(5) + t152 * t165;
t140 = -qJD(2) * pkin(2) + t161;
t130 = (pkin(5) * t156 + pkin(9) * t153) * qJD(2) + t132;
t129 = qJD(5) * pkin(9) + t167;
t128 = -qJD(5) * pkin(5) - t160;
t1 = [t168, t169, t157 * t162, -t154 * t162, -t140 * qJD(2), t143 * qJD(2), t143 ^ 2 / 0.2e1 + t150 ^ 2 * t168 + t140 ^ 2 / 0.2e1, -t134 * qJD(2), t135 * qJD(2), t135 ^ 2 / 0.2e1 + t134 ^ 2 / 0.2e1 + t145 ^ 2 / 0.2e1, t153 ^ 2 * t169, t153 * t158 * t156, -t153 * t163, -t156 * t163, qJD(5) ^ 2 / 0.2e1, qJD(5) * t160 + t132 * t164, -qJD(5) * t167 - t132 * t165, t142 ^ 2 / 0.2e1, t142 * t141, t142 * t146, t141 * t146, t146 ^ 2 / 0.2e1 (-t152 * t129 + t155 * t130) * t146 - t128 * t141 -(t155 * t129 + t152 * t130) * t146 + t128 * t142;];
T_reg  = t1;
