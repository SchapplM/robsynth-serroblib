% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:48:27
% EndTime: 2019-03-09 04:48:27
% DurationCPUTime: 0.08s
% Computational Cost: add. (324->43), mult. (615->88), div. (0->0), fcn. (357->6), ass. (0->33)
t160 = qJD(1) ^ 2;
t168 = t160 / 0.2e1;
t167 = cos(qJ(4));
t166 = t160 * qJ(2);
t157 = sin(qJ(4));
t159 = cos(qJ(3));
t164 = qJD(1) * t159;
t150 = t157 * qJD(3) + t167 * t164;
t158 = sin(qJ(3));
t152 = t158 * qJD(1) + qJD(4);
t146 = (pkin(3) * t158 - pkin(8) * t159 + qJ(2)) * qJD(1);
t151 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t147 = qJD(3) * pkin(8) + t158 * t151;
t161 = t167 * t146 - t157 * t147;
t136 = t152 * pkin(4) - t150 * qJ(5) + t161;
t149 = -t167 * qJD(3) + t157 * t164;
t165 = t157 * t146 + t167 * t147;
t138 = -t149 * qJ(5) + t165;
t155 = sin(pkin(9));
t156 = cos(pkin(9));
t133 = t155 * t136 + t156 * t138;
t163 = qJD(3) * t151;
t162 = qJD(1) * qJD(3);
t148 = -qJD(3) * pkin(3) - t159 * t151;
t132 = t156 * t136 - t155 * t138;
t141 = t149 * pkin(4) + qJD(5) + t148;
t153 = -qJD(1) * pkin(1) + qJD(2);
t140 = -t155 * t149 + t156 * t150;
t139 = t156 * t149 + t155 * t150;
t134 = t139 * pkin(5) - t140 * qJ(6) + t141;
t131 = t152 * qJ(6) + t133;
t130 = -t152 * pkin(5) + qJD(6) - t132;
t1 = [t168, 0, 0, t153 * qJD(1), t166, qJ(2) ^ 2 * t168 + t153 ^ 2 / 0.2e1, t159 ^ 2 * t168, -t159 * t160 * t158, t159 * t162, -t158 * t162, qJD(3) ^ 2 / 0.2e1, t158 * t166 + t159 * t163, -t158 * t163 + t159 * t166, t150 ^ 2 / 0.2e1, -t150 * t149, t150 * t152, -t149 * t152, t152 ^ 2 / 0.2e1, t148 * t149 + t161 * t152, t148 * t150 - t165 * t152, -t132 * t140 - t133 * t139, t133 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1 + t141 ^ 2 / 0.2e1, -t130 * t152 + t134 * t139, t130 * t140 - t131 * t139, t131 * t152 - t134 * t140, t131 ^ 2 / 0.2e1 + t134 ^ 2 / 0.2e1 + t130 ^ 2 / 0.2e1;];
T_reg  = t1;
