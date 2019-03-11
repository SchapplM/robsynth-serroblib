% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:44:57
% EndTime: 2019-03-09 01:44:57
% DurationCPUTime: 0.08s
% Computational Cost: add. (175->40), mult. (357->85), div. (0->0), fcn. (197->8), ass. (0->34)
t143 = sin(pkin(10));
t145 = cos(pkin(10));
t148 = sin(qJ(4));
t150 = cos(qJ(4));
t133 = (t143 * t150 + t145 * t148) * qJD(1);
t144 = sin(pkin(9));
t137 = (pkin(1) * t144 + qJ(3)) * qJD(1);
t151 = qJD(1) ^ 2;
t161 = t151 / 0.2e1;
t146 = cos(pkin(9));
t156 = -pkin(1) * t146 - pkin(2);
t135 = qJD(3) + (-pkin(7) + t156) * qJD(1);
t154 = -t148 * qJD(2) + t150 * t135;
t125 = -t150 * qJD(1) * qJ(5) + qJD(4) * pkin(4) + t154;
t159 = qJD(1) * t148;
t160 = t150 * qJD(2) + t148 * t135;
t126 = -qJ(5) * t159 + t160;
t121 = t143 * t125 + t145 * t126;
t158 = t137 * qJD(1);
t157 = qJD(1) * qJD(4);
t120 = t145 * t125 - t143 * t126;
t132 = pkin(4) * t159 + qJD(5) + t137;
t149 = cos(qJ(6));
t147 = sin(qJ(6));
t142 = qJD(2) ^ 2 / 0.2e1;
t136 = t156 * qJD(1) + qJD(3);
t134 = (-t143 * t148 + t145 * t150) * qJD(1);
t129 = qJD(6) + t133;
t128 = t147 * qJD(4) + t149 * t134;
t127 = -t149 * qJD(4) + t147 * t134;
t122 = t133 * pkin(5) - t134 * pkin(8) + t132;
t119 = qJD(4) * pkin(8) + t121;
t118 = -qJD(4) * pkin(5) - t120;
t1 = [t161, 0, 0, t142 + (t144 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t151, t136 * qJD(1), t158, t142 + t137 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1, t150 ^ 2 * t161, -t150 * t151 * t148, t150 * t157, -t148 * t157, qJD(4) ^ 2 / 0.2e1, t154 * qJD(4) + t148 * t158, -t160 * qJD(4) + t150 * t158, -t120 * t134 - t121 * t133, t121 ^ 2 / 0.2e1 + t120 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1, t128 ^ 2 / 0.2e1, -t128 * t127, t128 * t129, -t127 * t129, t129 ^ 2 / 0.2e1 (-t147 * t119 + t149 * t122) * t129 + t118 * t127 -(t149 * t119 + t147 * t122) * t129 + t118 * t128;];
T_reg  = t1;
