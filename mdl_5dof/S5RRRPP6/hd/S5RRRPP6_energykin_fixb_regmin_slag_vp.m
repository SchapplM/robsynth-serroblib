% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:01:59
% EndTime: 2019-12-31 21:01:59
% DurationCPUTime: 0.08s
% Computational Cost: add. (259->37), mult. (574->80), div. (0->0), fcn. (357->6), ass. (0->33)
t150 = qJD(1) ^ 2;
t160 = t150 / 0.2e1;
t159 = cos(qJ(3));
t149 = cos(qJ(2));
t158 = t149 * t150;
t147 = sin(qJ(3));
t148 = sin(qJ(2));
t156 = qJD(1) * t148;
t137 = t147 * qJD(2) + t159 * t156;
t155 = t149 * qJD(1);
t141 = -qJD(3) + t155;
t135 = (-pkin(2) * t149 - pkin(7) * t148 - pkin(1)) * qJD(1);
t140 = pkin(6) * t155 + qJD(2) * pkin(7);
t151 = t159 * t135 - t147 * t140;
t126 = -t141 * pkin(3) - t137 * qJ(4) + t151;
t136 = -t159 * qJD(2) + t147 * t156;
t157 = t147 * t135 + t159 * t140;
t128 = -t136 * qJ(4) + t157;
t145 = sin(pkin(8));
t146 = cos(pkin(8));
t123 = t145 * t126 + t146 * t128;
t154 = qJD(1) * qJD(2);
t153 = t148 * t154;
t152 = t149 * t154;
t139 = -qJD(2) * pkin(2) + pkin(6) * t156;
t122 = t146 * t126 - t145 * t128;
t131 = t136 * pkin(3) + qJD(4) + t139;
t130 = -t145 * t136 + t146 * t137;
t129 = t146 * t136 + t145 * t137;
t124 = t129 * pkin(4) - t130 * qJ(5) + t131;
t121 = -t141 * qJ(5) + t123;
t120 = t141 * pkin(4) + qJD(5) - t122;
t1 = [t160, 0, 0, t148 ^ 2 * t160, t148 * t158, t153, t152, qJD(2) ^ 2 / 0.2e1, pkin(1) * t158 - pkin(6) * t153, -t150 * pkin(1) * t148 - pkin(6) * t152, t137 ^ 2 / 0.2e1, -t137 * t136, -t137 * t141, t136 * t141, t141 ^ 2 / 0.2e1, t139 * t136 - t151 * t141, t139 * t137 + t157 * t141, -t122 * t130 - t123 * t129, t123 ^ 2 / 0.2e1 + t122 ^ 2 / 0.2e1 + t131 ^ 2 / 0.2e1, t120 * t141 + t124 * t129, t120 * t130 - t121 * t129, -t121 * t141 - t124 * t130, t121 ^ 2 / 0.2e1 + t124 ^ 2 / 0.2e1 + t120 ^ 2 / 0.2e1;];
T_reg = t1;
