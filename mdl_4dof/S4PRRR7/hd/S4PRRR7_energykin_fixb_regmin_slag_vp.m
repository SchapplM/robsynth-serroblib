% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:47
% EndTime: 2019-12-31 16:36:47
% DurationCPUTime: 0.05s
% Computational Cost: add. (43->21), mult. (130->56), div. (0->0), fcn. (83->8), ass. (0->26)
t109 = qJD(2) ^ 2;
t119 = t109 / 0.2e1;
t104 = sin(qJ(3));
t107 = cos(qJ(3));
t116 = qJD(1) * cos(pkin(4));
t105 = sin(qJ(2));
t117 = qJD(1) * sin(pkin(4));
t95 = qJD(2) * pkin(6) + t105 * t117;
t118 = t104 * t116 + t107 * t95;
t115 = qJD(2) * t104;
t114 = t107 * qJD(2);
t113 = qJD(2) * qJD(3);
t108 = cos(qJ(2));
t112 = t108 * t117;
t111 = qJD(2) * t117;
t110 = -t104 * t95 + t107 * t116;
t106 = cos(qJ(4));
t103 = sin(qJ(4));
t98 = -qJD(4) + t114;
t96 = -qJD(2) * pkin(2) - t112;
t94 = t103 * qJD(3) + t106 * t115;
t93 = -t106 * qJD(3) + t103 * t115;
t91 = -t112 + (-pkin(3) * t107 - pkin(7) * t104 - pkin(2)) * qJD(2);
t90 = qJD(3) * pkin(7) + t118;
t89 = -qJD(3) * pkin(3) - t110;
t1 = [qJD(1) ^ 2 / 0.2e1, t119, t108 * t111, -t105 * t111, t104 ^ 2 * t119, t104 * t109 * t107, t104 * t113, t107 * t113, qJD(3) ^ 2 / 0.2e1, t110 * qJD(3) - t96 * t114, -t118 * qJD(3) + t96 * t115, t94 ^ 2 / 0.2e1, -t94 * t93, -t94 * t98, t93 * t98, t98 ^ 2 / 0.2e1, -(-t103 * t90 + t106 * t91) * t98 + t89 * t93, (t103 * t91 + t106 * t90) * t98 + t89 * t94;];
T_reg = t1;
