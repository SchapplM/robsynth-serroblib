% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:39
% EndTime: 2019-12-05 17:16:39
% DurationCPUTime: 0.09s
% Computational Cost: add. (130->35), mult. (321->81), div. (0->0), fcn. (231->10), ass. (0->36)
t153 = qJD(2) ^ 2;
t164 = t153 / 0.2e1;
t148 = sin(qJ(2));
t161 = qJD(1) * sin(pkin(5));
t136 = qJD(2) * pkin(7) + t148 * t161;
t151 = cos(qJ(3));
t160 = qJD(1) * cos(pkin(5));
t139 = t151 * t160;
t147 = sin(qJ(3));
t127 = qJD(3) * pkin(3) + t139 + (-pkin(8) * qJD(2) - t136) * t147;
t158 = qJD(2) * t151;
t162 = t151 * t136 + t147 * t160;
t128 = pkin(8) * t158 + t162;
t146 = sin(qJ(4));
t150 = cos(qJ(4));
t163 = t146 * t127 + t150 * t128;
t159 = qJD(2) * t147;
t157 = qJD(2) * qJD(3);
t152 = cos(qJ(2));
t156 = t152 * t161;
t155 = qJD(2) * t161;
t133 = t146 * t159 - t150 * t158;
t154 = t150 * t127 - t146 * t128;
t131 = -t156 + (-pkin(3) * t151 - pkin(2)) * qJD(2);
t149 = cos(qJ(5));
t145 = sin(qJ(5));
t142 = qJD(3) + qJD(4);
t137 = -qJD(2) * pkin(2) - t156;
t134 = (t146 * t151 + t147 * t150) * qJD(2);
t132 = qJD(5) + t133;
t130 = t149 * t134 + t145 * t142;
t129 = t145 * t134 - t149 * t142;
t124 = t133 * pkin(4) - t134 * pkin(9) + t131;
t123 = t142 * pkin(9) + t163;
t122 = -t142 * pkin(4) - t154;
t1 = [qJD(1) ^ 2 / 0.2e1, t164, t152 * t155, -t148 * t155, t147 ^ 2 * t164, t147 * t153 * t151, t147 * t157, t151 * t157, qJD(3) ^ 2 / 0.2e1, (-t147 * t136 + t139) * qJD(3) - t137 * t158, -t162 * qJD(3) + t137 * t159, t134 ^ 2 / 0.2e1, -t134 * t133, t134 * t142, -t133 * t142, t142 ^ 2 / 0.2e1, t131 * t133 + t154 * t142, t131 * t134 - t163 * t142, t130 ^ 2 / 0.2e1, -t130 * t129, t130 * t132, -t129 * t132, t132 ^ 2 / 0.2e1, (-t145 * t123 + t149 * t124) * t132 + t122 * t129, -(t149 * t123 + t145 * t124) * t132 + t122 * t130;];
T_reg = t1;
