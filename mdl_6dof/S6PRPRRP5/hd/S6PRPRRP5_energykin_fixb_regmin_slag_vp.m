% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:16:42
% EndTime: 2019-03-08 20:16:42
% DurationCPUTime: 0.07s
% Computational Cost: add. (127->36), mult. (272->76), div. (0->0), fcn. (164->8), ass. (0->36)
t159 = qJD(2) ^ 2;
t176 = t159 / 0.2e1;
t175 = qJD(1) ^ 2 / 0.2e1;
t174 = cos(qJ(5));
t158 = cos(qJ(2));
t171 = qJD(1) * sin(pkin(6));
t161 = -t158 * t171 + qJD(3);
t142 = (-pkin(2) - pkin(8)) * qJD(2) + t161;
t155 = sin(qJ(4));
t157 = cos(qJ(4));
t153 = cos(pkin(6));
t170 = qJD(1) * t153;
t172 = t155 * t142 + t157 * t170;
t137 = qJD(4) * pkin(9) + t172;
t156 = sin(qJ(2));
t165 = t156 * t171;
t140 = t165 + (pkin(4) * t155 - pkin(9) * t157 + qJ(3)) * qJD(2);
t154 = sin(qJ(5));
t173 = t174 * t137 + t154 * t140;
t169 = qJD(2) * t157;
t146 = qJD(2) * qJ(3) + t165;
t168 = t146 * qJD(2);
t167 = t155 * qJD(2);
t166 = qJD(2) * qJD(4);
t164 = qJD(2) * t171;
t163 = -t154 * t137 + t174 * t140;
t162 = t157 * t142 - t155 * t170;
t136 = -qJD(4) * pkin(4) - t162;
t150 = qJD(5) + t167;
t145 = t154 * qJD(4) + t174 * t169;
t144 = -t174 * qJD(4) + t154 * t169;
t143 = -qJD(2) * pkin(2) + t161;
t134 = t144 * pkin(5) + qJD(6) + t136;
t133 = -t144 * qJ(6) + t173;
t132 = t150 * pkin(5) - t145 * qJ(6) + t163;
t1 = [t175, t176, t158 * t164, -t156 * t164, t143 * qJD(2), t168, t153 ^ 2 * t175 + t146 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1, t157 ^ 2 * t176, -t157 * t159 * t155, t157 * t166, -t155 * t166, qJD(4) ^ 2 / 0.2e1, t162 * qJD(4) + t146 * t167, -t172 * qJD(4) + t157 * t168, t145 ^ 2 / 0.2e1, -t145 * t144, t145 * t150, -t144 * t150, t150 ^ 2 / 0.2e1, t136 * t144 + t163 * t150, t136 * t145 - t173 * t150, -t132 * t145 - t133 * t144, t133 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1 + t134 ^ 2 / 0.2e1;];
T_reg  = t1;
