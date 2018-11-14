% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4RPPP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:35
% EndTime: 2018-11-14 13:45:35
% DurationCPUTime: 0.10s
% Computational Cost: add. (87->34), mult. (287->72), div. (0->0), fcn. (168->4), ass. (0->24)
t101 = sin(pkin(6));
t103 = cos(pkin(6));
t102 = sin(pkin(4));
t110 = qJD(1) * t102;
t106 = qJ(2) * t110;
t104 = cos(pkin(4));
t109 = qJD(1) * t104;
t108 = pkin(1) * t109;
t96 = t101 * t108 + t103 * t106;
t113 = t101 * t102;
t112 = t102 * t103;
t97 = t101 * t106;
t111 = qJD(3) + t97;
t107 = -pkin(1) * t103 - pkin(2);
t105 = -qJ(3) * t101 - pkin(1);
t100 = -pkin(1) * t110 + qJD(2);
t95 = t103 * t108 - t97;
t94 = -qJ(3) * t109 - t96;
t93 = qJD(2) + (-pkin(2) * t103 + t105) * t110;
t92 = t107 * t109 + t111;
t91 = qJD(2) + ((-pkin(2) - qJ(4)) * t103 + t105) * t110;
t90 = qJD(4) + (pkin(3) * t112 + qJ(3) * t104) * qJD(1) + t96;
t89 = (pkin(3) * t113 + (-qJ(4) + t107) * t104) * qJD(1) + t111;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0 (-t100 * t112 + t104 * t95) * qJD(1) (t100 * t113 - t104 * t96) * qJD(1) (-t101 * t95 + t103 * t96) * t110, t96 ^ 2 / 0.2e1 + t95 ^ 2 / 0.2e1 + t100 ^ 2 / 0.2e1 (t101 * t92 - t103 * t94) * t110 (t104 * t92 + t93 * t112) * qJD(1) (-t104 * t94 - t93 * t113) * qJD(1), t93 ^ 2 / 0.2e1 + t94 ^ 2 / 0.2e1 + t92 ^ 2 / 0.2e1 (t101 * t89 + t103 * t90) * t110 (t104 * t90 - t91 * t113) * qJD(1) (-t104 * t89 - t91 * t112) * qJD(1), t91 ^ 2 / 0.2e1 + t89 ^ 2 / 0.2e1 + t90 ^ 2 / 0.2e1;];
T_reg  = t1;
