% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% T_reg [1x31]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:24:33
% EndTime: 2019-03-09 06:24:33
% DurationCPUTime: 0.12s
% Computational Cost: add. (320->43), mult. (602->90), div. (0->0), fcn. (371->6), ass. (0->34)
t149 = qJD(1) ^ 2;
t158 = t149 / 0.2e1;
t157 = t149 * qJ(2);
t141 = qJD(3) + qJD(4);
t148 = cos(qJ(3));
t137 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t152 = -pkin(8) * qJD(1) + t137;
t131 = qJD(3) * pkin(3) + t152 * t148;
t145 = sin(qJ(3));
t132 = t152 * t145;
t144 = sin(qJ(4));
t147 = cos(qJ(4));
t155 = t144 * t131 + t147 * t132;
t123 = t141 * pkin(9) + t155;
t134 = (t144 * t148 + t145 * t147) * qJD(1);
t135 = (-t144 * t145 + t147 * t148) * qJD(1);
t136 = (pkin(3) * t145 + qJ(2)) * qJD(1);
t125 = t134 * pkin(4) - t135 * pkin(9) + t136;
t143 = sin(qJ(5));
t146 = cos(qJ(5));
t156 = t146 * t123 + t143 * t125;
t154 = qJD(3) * t137;
t153 = qJD(1) * qJD(3);
t151 = t147 * t131 - t144 * t132;
t150 = -t143 * t123 + t146 * t125;
t122 = -t141 * pkin(4) - t151;
t140 = -qJD(1) * pkin(1) + qJD(2);
t133 = qJD(5) + t134;
t127 = t146 * t135 + t143 * t141;
t126 = t143 * t135 - t146 * t141;
t120 = t126 * pkin(5) - t127 * qJ(6) + t122;
t119 = t133 * qJ(6) + t156;
t118 = -t133 * pkin(5) + qJD(6) - t150;
t1 = [t158, 0, 0, t140 * qJD(1), t157, qJ(2) ^ 2 * t158 + t140 ^ 2 / 0.2e1, t148 ^ 2 * t158, -t148 * t149 * t145, t148 * t153, -t145 * t153, qJD(3) ^ 2 / 0.2e1, t145 * t157 + t148 * t154, -t145 * t154 + t148 * t157, t135 ^ 2 / 0.2e1, -t135 * t134, t135 * t141, -t134 * t141, t141 ^ 2 / 0.2e1, t136 * t134 + t151 * t141, t136 * t135 - t155 * t141, t127 ^ 2 / 0.2e1, -t127 * t126, t127 * t133, -t126 * t133, t133 ^ 2 / 0.2e1, t122 * t126 + t150 * t133, t122 * t127 - t156 * t133, -t118 * t133 + t120 * t126, t118 * t127 - t119 * t126, t119 * t133 - t120 * t127, t119 ^ 2 / 0.2e1 + t120 ^ 2 / 0.2e1 + t118 ^ 2 / 0.2e1;];
T_reg  = t1;
