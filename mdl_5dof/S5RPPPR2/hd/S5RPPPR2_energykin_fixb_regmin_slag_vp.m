% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:32
% EndTime: 2022-01-23 08:59:32
% DurationCPUTime: 0.09s
% Computational Cost: add. (208->43), mult. (553->95), div. (0->0), fcn. (353->8), ass. (0->35)
t159 = sin(pkin(7));
t161 = cos(pkin(8));
t171 = t159 * t161;
t162 = cos(pkin(7));
t147 = qJD(2) + (-pkin(2) * t162 - qJ(3) * t159 - pkin(1)) * qJD(1);
t158 = sin(pkin(8));
t169 = qJD(1) * t162;
t168 = qJ(2) * t169;
t140 = t158 * t147 + t161 * t168;
t136 = -qJ(4) * t169 + t140;
t170 = qJD(1) * t159;
t151 = qJ(2) * t170 + qJD(3);
t142 = (pkin(3) * t158 - qJ(4) * t161) * t170 + t151;
t157 = sin(pkin(9));
t160 = cos(pkin(9));
t133 = t160 * t136 + t157 * t142;
t167 = t158 * t170;
t139 = t161 * t147 - t158 * t168;
t132 = -t157 * t136 + t160 * t142;
t135 = pkin(3) * t169 + qJD(4) - t139;
t145 = (t157 * t171 + t160 * t162) * qJD(1);
t165 = qJD(1) ^ 2;
t164 = cos(qJ(5));
t163 = sin(qJ(5));
t156 = t162 ^ 2;
t155 = t159 ^ 2;
t154 = -qJD(1) * pkin(1) + qJD(2);
t146 = (-t157 * t162 + t160 * t171) * qJD(1);
t143 = qJD(5) + t145;
t138 = t164 * t146 + t163 * t167;
t137 = t163 * t146 - t164 * t167;
t131 = pkin(6) * t167 + t133;
t130 = -pkin(4) * t167 - t132;
t129 = t145 * pkin(4) - t146 * pkin(6) + t135;
t1 = [t165 / 0.2e1, 0, 0, -t154 * t169, t154 * t170, (t155 + t156) * t165 * qJ(2), t154 ^ 2 / 0.2e1 + (t156 / 0.2e1 + t155 / 0.2e1) * qJ(2) ^ 2 * t165, (t151 * t158 * t159 - t139 * t162) * qJD(1), (t140 * t162 + t151 * t171) * qJD(1), (-t139 * t161 - t140 * t158) * t170, t140 ^ 2 / 0.2e1 + t139 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1, t132 * t167 + t135 * t145, -t133 * t167 + t135 * t146, -t132 * t146 - t133 * t145, t133 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1 + t135 ^ 2 / 0.2e1, t138 ^ 2 / 0.2e1, -t138 * t137, t138 * t143, -t137 * t143, t143 ^ 2 / 0.2e1, (t164 * t129 - t163 * t131) * t143 + t130 * t137, -(t163 * t129 + t164 * t131) * t143 + t130 * t138;];
T_reg = t1;
