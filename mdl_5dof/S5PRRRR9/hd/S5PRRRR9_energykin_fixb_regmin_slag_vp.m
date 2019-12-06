% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRR9
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:21:10
% EndTime: 2019-12-05 17:21:10
% DurationCPUTime: 0.08s
% Computational Cost: add. (128->34), mult. (312->79), div. (0->0), fcn. (222->10), ass. (0->36)
t158 = qJD(2) ^ 2;
t171 = t158 / 0.2e1;
t170 = cos(qJ(4));
t154 = sin(qJ(2));
t167 = qJD(1) * sin(pkin(5));
t141 = qJD(2) * pkin(7) + t154 * t167;
t153 = sin(qJ(3));
t156 = cos(qJ(3));
t166 = qJD(1) * cos(pkin(5));
t168 = t156 * t141 + t153 * t166;
t132 = qJD(3) * pkin(8) + t168;
t157 = cos(qJ(2));
t162 = t157 * t167;
t135 = -t162 + (-pkin(3) * t156 - pkin(8) * t153 - pkin(2)) * qJD(2);
t152 = sin(qJ(4));
t169 = t170 * t132 + t152 * t135;
t165 = qJD(2) * t153;
t164 = t156 * qJD(2);
t163 = qJD(2) * qJD(3);
t161 = qJD(2) * t167;
t160 = -t152 * t132 + t170 * t135;
t146 = -qJD(4) + t164;
t159 = -t153 * t141 + t156 * t166;
t131 = -qJD(3) * pkin(3) - t159;
t155 = cos(qJ(5));
t151 = sin(qJ(5));
t143 = -qJD(5) + t146;
t142 = -qJD(2) * pkin(2) - t162;
t140 = t152 * qJD(3) + t170 * t165;
t139 = -t170 * qJD(3) + t152 * t165;
t129 = -t151 * t139 + t155 * t140;
t128 = t155 * t139 + t151 * t140;
t127 = t139 * pkin(4) + t131;
t126 = -t139 * pkin(9) + t169;
t125 = -t146 * pkin(4) - t140 * pkin(9) + t160;
t1 = [qJD(1) ^ 2 / 0.2e1, t171, t157 * t161, -t154 * t161, t153 ^ 2 * t171, t153 * t158 * t156, t153 * t163, t156 * t163, qJD(3) ^ 2 / 0.2e1, t159 * qJD(3) - t142 * t164, -t168 * qJD(3) + t142 * t165, t140 ^ 2 / 0.2e1, -t140 * t139, -t140 * t146, t139 * t146, t146 ^ 2 / 0.2e1, t131 * t139 - t160 * t146, t131 * t140 + t169 * t146, t129 ^ 2 / 0.2e1, -t129 * t128, -t129 * t143, t128 * t143, t143 ^ 2 / 0.2e1, t127 * t128 - (t155 * t125 - t151 * t126) * t143, t127 * t129 + (t151 * t125 + t155 * t126) * t143;];
T_reg = t1;
