% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:44
% EndTime: 2019-12-31 19:10:44
% DurationCPUTime: 0.13s
% Computational Cost: add. (212->42), mult. (554->88), div. (0->0), fcn. (402->8), ass. (0->37)
t162 = cos(qJ(4));
t161 = pkin(6) + qJ(2);
t150 = sin(qJ(3));
t152 = cos(qJ(3));
t147 = cos(pkin(9));
t157 = qJD(1) * t147;
t146 = sin(pkin(9));
t158 = qJD(1) * t146;
t136 = t150 * t158 - t152 * t157;
t137 = (t146 * t152 + t147 * t150) * qJD(1);
t140 = qJD(2) + (-pkin(2) * t147 - pkin(1)) * qJD(1);
t124 = t136 * pkin(3) - t137 * pkin(7) + t140;
t138 = t161 * t158;
t139 = t161 * t157;
t159 = -t150 * t138 + t152 * t139;
t127 = qJD(3) * pkin(7) + t159;
t149 = sin(qJ(4));
t160 = t149 * t124 + t162 * t127;
t156 = t162 * t124 - t149 * t127;
t155 = -t152 * t138 - t150 * t139;
t132 = qJD(4) + t136;
t126 = -qJD(3) * pkin(3) - t155;
t153 = qJD(1) ^ 2;
t151 = cos(qJ(5));
t148 = sin(qJ(5));
t145 = t147 ^ 2;
t144 = t146 ^ 2;
t142 = -qJD(1) * pkin(1) + qJD(2);
t131 = qJD(5) + t132;
t130 = t149 * qJD(3) + t162 * t137;
t129 = -t162 * qJD(3) + t149 * t137;
t121 = -t148 * t129 + t151 * t130;
t120 = t151 * t129 + t148 * t130;
t119 = t129 * pkin(4) + t126;
t118 = -t129 * pkin(8) + t160;
t117 = t132 * pkin(4) - t130 * pkin(8) + t156;
t1 = [t153 / 0.2e1, 0, 0, -t142 * t157, t142 * t158, (t144 + t145) * t153 * qJ(2), t142 ^ 2 / 0.2e1 + (t145 / 0.2e1 + t144 / 0.2e1) * qJ(2) ^ 2 * t153, t137 ^ 2 / 0.2e1, -t137 * t136, t137 * qJD(3), -t136 * qJD(3), qJD(3) ^ 2 / 0.2e1, t155 * qJD(3) + t140 * t136, -t159 * qJD(3) + t140 * t137, t130 ^ 2 / 0.2e1, -t130 * t129, t130 * t132, -t129 * t132, t132 ^ 2 / 0.2e1, t126 * t129 + t156 * t132, t126 * t130 - t160 * t132, t121 ^ 2 / 0.2e1, -t121 * t120, t121 * t131, -t120 * t131, t131 ^ 2 / 0.2e1, (t151 * t117 - t148 * t118) * t131 + t119 * t120, -(t148 * t117 + t151 * t118) * t131 + t119 * t121;];
T_reg = t1;
