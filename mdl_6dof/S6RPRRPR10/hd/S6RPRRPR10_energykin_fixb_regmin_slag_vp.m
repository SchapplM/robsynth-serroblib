% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% 
% Output:
% T_reg [1x31]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:37:09
% EndTime: 2019-03-09 05:37:09
% DurationCPUTime: 0.13s
% Computational Cost: add. (225->45), mult. (427->91), div. (0->0), fcn. (241->6), ass. (0->35)
t169 = -pkin(4) - pkin(5);
t159 = qJD(1) ^ 2;
t168 = t159 / 0.2e1;
t167 = t159 * qJ(2);
t155 = sin(qJ(3));
t158 = cos(qJ(3));
t141 = (pkin(3) * t155 - pkin(8) * t158 + qJ(2)) * qJD(1);
t149 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t142 = qJD(3) * pkin(8) + t149 * t155;
t154 = sin(qJ(4));
t157 = cos(qJ(4));
t166 = t154 * t141 + t157 * t142;
t165 = qJD(1) * t158;
t164 = qJD(3) * t149;
t163 = qJD(1) * qJD(3);
t150 = qJD(1) * t155 + qJD(4);
t133 = t150 * qJ(5) + t166;
t162 = t157 * t141 - t154 * t142;
t161 = qJD(5) - t162;
t143 = -qJD(3) * pkin(3) - t158 * t149;
t145 = qJD(3) * t154 + t157 * t165;
t160 = t145 * qJ(5) - t143;
t156 = cos(qJ(6));
t153 = sin(qJ(6));
t151 = -qJD(1) * pkin(1) + qJD(2);
t147 = -qJD(6) + t150;
t144 = -t157 * qJD(3) + t154 * t165;
t136 = t144 * t153 + t145 * t156;
t135 = -t156 * t144 + t145 * t153;
t134 = pkin(4) * t144 - t160;
t132 = -pkin(4) * t150 + t161;
t131 = t144 * t169 + t160;
t130 = pkin(9) * t144 + t133;
t129 = -t145 * pkin(9) + t150 * t169 + t161;
t1 = [t168, 0, 0, t151 * qJD(1), t167, qJ(2) ^ 2 * t168 + t151 ^ 2 / 0.2e1, t158 ^ 2 * t168, -t158 * t159 * t155, t158 * t163, -t155 * t163, qJD(3) ^ 2 / 0.2e1, t155 * t167 + t158 * t164, -t155 * t164 + t158 * t167, t145 ^ 2 / 0.2e1, -t145 * t144, t145 * t150, -t144 * t150, t150 ^ 2 / 0.2e1, t143 * t144 + t150 * t162, t143 * t145 - t150 * t166, -t132 * t150 + t134 * t144, t132 * t145 - t133 * t144, t133 * t150 - t134 * t145, t133 ^ 2 / 0.2e1 + t134 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1, t136 ^ 2 / 0.2e1, -t136 * t135, -t136 * t147, t135 * t147, t147 ^ 2 / 0.2e1 -(t129 * t156 - t130 * t153) * t147 + t131 * t135 (t129 * t153 + t130 * t156) * t147 + t131 * t136;];
T_reg  = t1;
