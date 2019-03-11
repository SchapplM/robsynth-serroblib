% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:57:25
% EndTime: 2019-03-09 02:57:25
% DurationCPUTime: 0.08s
% Computational Cost: add. (226->43), mult. (453->88), div. (0->0), fcn. (253->6), ass. (0->35)
t163 = pkin(4) + pkin(8);
t154 = qJD(1) ^ 2;
t162 = t154 / 0.2e1;
t161 = t154 * qJ(2);
t153 = cos(qJ(3));
t143 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t157 = -qJ(4) * qJD(1) + t143;
t136 = qJD(3) * pkin(3) + t157 * t153;
t151 = sin(qJ(3));
t138 = t157 * t151;
t148 = sin(pkin(9));
t149 = cos(pkin(9));
t129 = t148 * t136 + t149 * t138;
t160 = qJD(1) * t151;
t159 = qJD(3) * t143;
t158 = qJD(1) * qJD(3);
t142 = pkin(3) * t160 + qJD(1) * qJ(2) + qJD(4);
t128 = t149 * t136 - t148 * t138;
t141 = t149 * t153 * qJD(1) - t148 * t160;
t127 = -qJD(3) * qJ(5) - t129;
t156 = qJD(5) - t128;
t155 = -t141 * qJ(5) + t142;
t152 = cos(qJ(6));
t150 = sin(qJ(6));
t146 = -qJD(1) * pkin(1) + qJD(2);
t140 = (t148 * t153 + t149 * t151) * qJD(1);
t139 = qJD(6) + t141;
t132 = t152 * qJD(3) + t150 * t140;
t131 = t150 * qJD(3) - t152 * t140;
t130 = t140 * pkin(4) + t155;
t126 = -qJD(3) * pkin(4) + t156;
t125 = t163 * t140 + t155;
t124 = -t140 * pkin(5) - t127;
t123 = t141 * pkin(5) - t163 * qJD(3) + t156;
t1 = [t162, 0, 0, t146 * qJD(1), t161, qJ(2) ^ 2 * t162 + t146 ^ 2 / 0.2e1, t153 ^ 2 * t162, -t153 * t154 * t151, t153 * t158, -t151 * t158, qJD(3) ^ 2 / 0.2e1, t151 * t161 + t153 * t159, -t151 * t159 + t153 * t161, -t128 * t141 - t129 * t140, t129 ^ 2 / 0.2e1 + t128 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1, t126 * t141 + t127 * t140, t126 * qJD(3) - t130 * t140, -t127 * qJD(3) - t130 * t141, t130 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1 + t126 ^ 2 / 0.2e1, t132 ^ 2 / 0.2e1, -t132 * t131, t132 * t139, -t131 * t139, t139 ^ 2 / 0.2e1 (t152 * t123 - t150 * t125) * t139 + t124 * t131 -(t150 * t123 + t152 * t125) * t139 + t124 * t132;];
T_reg  = t1;
