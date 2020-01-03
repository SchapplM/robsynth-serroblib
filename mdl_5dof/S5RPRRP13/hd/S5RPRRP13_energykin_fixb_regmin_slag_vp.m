% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP13_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:36
% EndTime: 2019-12-31 18:59:36
% DurationCPUTime: 0.06s
% Computational Cost: add. (126->31), mult. (247->68), div. (0->0), fcn. (118->4), ass. (0->24)
t105 = qJD(1) ^ 2;
t112 = t105 / 0.2e1;
t101 = sin(qJ(4));
t103 = cos(qJ(4));
t102 = sin(qJ(3));
t104 = cos(qJ(3));
t92 = (pkin(3) * t102 - pkin(7) * t104 + qJ(2)) * qJD(1);
t97 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t93 = qJD(3) * pkin(7) + t102 * t97;
t111 = t101 * t92 + t103 * t93;
t110 = qJD(3) * t97;
t109 = t105 * qJ(2);
t108 = qJD(1) * t104;
t107 = qJD(1) * qJD(3);
t94 = -qJD(3) * pkin(3) - t104 * t97;
t106 = -t101 * t93 + t103 * t92;
t99 = -qJD(1) * pkin(1) + qJD(2);
t98 = t102 * qJD(1) + qJD(4);
t96 = t101 * qJD(3) + t103 * t108;
t95 = -t103 * qJD(3) + t101 * t108;
t89 = t95 * pkin(4) - t96 * qJ(5) + t94;
t88 = t98 * qJ(5) + t111;
t87 = -t98 * pkin(4) + qJD(5) - t106;
t1 = [t112, 0, 0, t99 * qJD(1), t109, qJ(2) ^ 2 * t112 + t99 ^ 2 / 0.2e1, t104 ^ 2 * t112, -t104 * t105 * t102, t104 * t107, -t102 * t107, qJD(3) ^ 2 / 0.2e1, t102 * t109 + t104 * t110, -t102 * t110 + t104 * t109, t96 ^ 2 / 0.2e1, -t96 * t95, t96 * t98, -t95 * t98, t98 ^ 2 / 0.2e1, t106 * t98 + t94 * t95, -t111 * t98 + t94 * t96, -t87 * t98 + t89 * t95, t87 * t96 - t88 * t95, t88 * t98 - t89 * t96, t88 ^ 2 / 0.2e1 + t89 ^ 2 / 0.2e1 + t87 ^ 2 / 0.2e1;];
T_reg = t1;
