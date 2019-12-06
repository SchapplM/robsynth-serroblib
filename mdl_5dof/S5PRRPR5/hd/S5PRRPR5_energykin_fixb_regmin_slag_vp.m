% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% T_reg [1x20]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:27:53
% EndTime: 2019-12-05 16:27:53
% DurationCPUTime: 0.08s
% Computational Cost: add. (127->34), mult. (331->78), div. (0->0), fcn. (236->10), ass. (0->35)
t162 = qJD(2) ^ 2;
t171 = t162 / 0.2e1;
t158 = sin(qJ(2));
t169 = qJD(1) * sin(pkin(5));
t146 = qJD(2) * pkin(7) + t158 * t169;
t160 = cos(qJ(3));
t168 = qJD(1) * cos(pkin(5));
t150 = t160 * t168;
t157 = sin(qJ(3));
t137 = qJD(3) * pkin(3) + t150 + (-qJ(4) * qJD(2) - t146) * t157;
t166 = qJD(2) * t160;
t170 = t160 * t146 + t157 * t168;
t138 = qJ(4) * t166 + t170;
t152 = sin(pkin(10));
t154 = cos(pkin(10));
t133 = t152 * t137 + t154 * t138;
t167 = qJD(2) * t157;
t165 = qJD(2) * qJD(3);
t161 = cos(qJ(2));
t164 = t161 * t169;
t163 = qJD(2) * t169;
t143 = -t152 * t167 + t154 * t166;
t132 = t154 * t137 - t152 * t138;
t141 = -t164 + qJD(4) + (-pkin(3) * t160 - pkin(2)) * qJD(2);
t159 = cos(qJ(5));
t156 = sin(qJ(5));
t147 = -qJD(2) * pkin(2) - t164;
t144 = (t152 * t160 + t154 * t157) * qJD(2);
t142 = qJD(5) - t143;
t140 = t156 * qJD(3) + t159 * t144;
t139 = -t159 * qJD(3) + t156 * t144;
t134 = -t143 * pkin(4) - t144 * pkin(8) + t141;
t131 = qJD(3) * pkin(8) + t133;
t130 = -qJD(3) * pkin(4) - t132;
t1 = [qJD(1) ^ 2 / 0.2e1, t171, t161 * t163, -t158 * t163, t157 ^ 2 * t171, t157 * t162 * t160, t157 * t165, t160 * t165, qJD(3) ^ 2 / 0.2e1, (-t157 * t146 + t150) * qJD(3) - t147 * t166, -t170 * qJD(3) + t147 * t167, -t132 * t144 + t133 * t143, t133 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1 + t141 ^ 2 / 0.2e1, t140 ^ 2 / 0.2e1, -t140 * t139, t140 * t142, -t139 * t142, t142 ^ 2 / 0.2e1, (-t156 * t131 + t159 * t134) * t142 + t130 * t139, -(t159 * t131 + t156 * t134) * t142 + t130 * t140;];
T_reg = t1;
