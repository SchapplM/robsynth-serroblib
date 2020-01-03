% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:31
% EndTime: 2019-12-31 20:01:31
% DurationCPUTime: 0.08s
% Computational Cost: add. (228->36), mult. (555->79), div. (0->0), fcn. (362->6), ass. (0->34)
t149 = qJD(1) ^ 2;
t159 = t149 / 0.2e1;
t158 = pkin(6) + qJ(3);
t148 = cos(qJ(2));
t157 = t148 * t149;
t143 = sin(pkin(8));
t144 = cos(pkin(8));
t154 = qJD(1) * t148;
t146 = sin(qJ(2));
t155 = qJD(1) * t146;
t134 = -t143 * t155 + t144 * t154;
t135 = (t143 * t148 + t144 * t146) * qJD(1);
t140 = qJD(3) + (-pkin(2) * t148 - pkin(1)) * qJD(1);
t124 = -t134 * pkin(3) - t135 * pkin(7) + t140;
t138 = qJD(2) * pkin(2) - t158 * t155;
t139 = t158 * t154;
t129 = t143 * t138 + t144 * t139;
t127 = qJD(2) * pkin(7) + t129;
t145 = sin(qJ(4));
t147 = cos(qJ(4));
t156 = t145 * t124 + t147 * t127;
t153 = qJD(1) * qJD(2);
t152 = t146 * t153;
t151 = t148 * t153;
t128 = t144 * t138 - t143 * t139;
t150 = t147 * t124 - t145 * t127;
t126 = -qJD(2) * pkin(3) - t128;
t133 = qJD(4) - t134;
t131 = t145 * qJD(2) + t147 * t135;
t130 = -t147 * qJD(2) + t145 * t135;
t122 = t130 * pkin(4) - t131 * qJ(5) + t126;
t121 = t133 * qJ(5) + t156;
t120 = -t133 * pkin(4) + qJD(5) - t150;
t1 = [t159, 0, 0, t146 ^ 2 * t159, t146 * t157, t152, t151, qJD(2) ^ 2 / 0.2e1, pkin(1) * t157 - pkin(6) * t152, -t149 * pkin(1) * t146 - pkin(6) * t151, -t128 * t135 + t129 * t134, t129 ^ 2 / 0.2e1 + t128 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1, t131 ^ 2 / 0.2e1, -t131 * t130, t131 * t133, -t130 * t133, t133 ^ 2 / 0.2e1, t126 * t130 + t150 * t133, t126 * t131 - t156 * t133, -t120 * t133 + t122 * t130, t120 * t131 - t121 * t130, t121 * t133 - t122 * t131, t121 ^ 2 / 0.2e1 + t122 ^ 2 / 0.2e1 + t120 ^ 2 / 0.2e1;];
T_reg = t1;
