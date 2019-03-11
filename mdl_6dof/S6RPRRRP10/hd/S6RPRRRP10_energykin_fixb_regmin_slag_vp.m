% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRP10
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
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:32:31
% EndTime: 2019-03-09 06:32:32
% DurationCPUTime: 0.13s
% Computational Cost: add. (319->44), mult. (603->91), div. (0->0), fcn. (368->6), ass. (0->34)
t154 = qJD(1) ^ 2;
t164 = t154 / 0.2e1;
t163 = cos(qJ(4));
t162 = t154 * qJ(2);
t150 = sin(qJ(4));
t153 = cos(qJ(3));
t159 = qJD(1) * t153;
t143 = t150 * qJD(3) + t163 * t159;
t151 = sin(qJ(3));
t146 = t151 * qJD(1) + qJD(4);
t139 = (pkin(3) * t151 - pkin(8) * t153 + qJ(2)) * qJD(1);
t145 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t140 = qJD(3) * pkin(8) + t151 * t145;
t156 = t163 * t139 - t150 * t140;
t129 = t146 * pkin(4) - t143 * pkin(9) + t156;
t142 = -t163 * qJD(3) + t150 * t159;
t160 = t150 * t139 + t163 * t140;
t131 = -t142 * pkin(9) + t160;
t149 = sin(qJ(5));
t152 = cos(qJ(5));
t161 = t149 * t129 + t152 * t131;
t158 = qJD(3) * t145;
t157 = qJD(1) * qJD(3);
t141 = -qJD(3) * pkin(3) - t153 * t145;
t155 = t152 * t129 - t149 * t131;
t134 = t142 * pkin(4) + t141;
t147 = -qJD(1) * pkin(1) + qJD(2);
t144 = qJD(5) + t146;
t133 = -t149 * t142 + t152 * t143;
t132 = t152 * t142 + t149 * t143;
t127 = t132 * pkin(5) - t133 * qJ(6) + t134;
t126 = t144 * qJ(6) + t161;
t125 = -t144 * pkin(5) + qJD(6) - t155;
t1 = [t164, 0, 0, t147 * qJD(1), t162, qJ(2) ^ 2 * t164 + t147 ^ 2 / 0.2e1, t153 ^ 2 * t164, -t153 * t154 * t151, t153 * t157, -t151 * t157, qJD(3) ^ 2 / 0.2e1, t151 * t162 + t153 * t158, -t151 * t158 + t153 * t162, t143 ^ 2 / 0.2e1, -t143 * t142, t143 * t146, -t142 * t146, t146 ^ 2 / 0.2e1, t141 * t142 + t156 * t146, t141 * t143 - t160 * t146, t133 ^ 2 / 0.2e1, -t133 * t132, t133 * t144, -t132 * t144, t144 ^ 2 / 0.2e1, t134 * t132 + t155 * t144, t134 * t133 - t161 * t144, -t125 * t144 + t127 * t132, t125 * t133 - t126 * t132, t126 * t144 - t127 * t133, t126 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1 + t125 ^ 2 / 0.2e1;];
T_reg  = t1;
