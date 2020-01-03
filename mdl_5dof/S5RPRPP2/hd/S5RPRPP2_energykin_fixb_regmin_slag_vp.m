% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:11
% EndTime: 2019-12-31 18:11:11
% DurationCPUTime: 0.07s
% Computational Cost: add. (95->29), mult. (218->71), div. (0->0), fcn. (88->4), ass. (0->25)
t94 = qJD(1) ^ 2;
t107 = t94 / 0.2e1;
t106 = pkin(3) + pkin(4);
t90 = sin(pkin(7));
t86 = (pkin(1) * t90 + pkin(6)) * qJD(1);
t92 = sin(qJ(3));
t93 = cos(qJ(3));
t105 = t92 * qJD(2) + t93 * t86;
t104 = qJD(1) * t92;
t103 = qJD(1) * t93;
t91 = cos(pkin(7));
t99 = -pkin(1) * t91 - pkin(2);
t102 = t99 * t94;
t101 = qJ(5) * qJD(1);
t100 = qJD(1) * qJD(3);
t82 = qJD(3) * qJ(4) + t105;
t98 = t93 * qJD(2) - t92 * t86;
t97 = qJD(4) - t98;
t96 = qJ(4) * t92 - t99;
t83 = (-pkin(3) * t93 - t96) * qJD(1);
t81 = -qJD(3) * pkin(3) + t97;
t80 = -t93 * t101 + t82;
t79 = qJD(5) + (t106 * t93 + t96) * qJD(1);
t78 = -t106 * qJD(3) - t92 * t101 + t97;
t1 = [t107, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t90 ^ 2 / 0.2e1 + t91 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t94, t92 ^ 2 * t107, t92 * t94 * t93, t92 * t100, t93 * t100, qJD(3) ^ 2 / 0.2e1, t98 * qJD(3) - t93 * t102, -t105 * qJD(3) + t92 * t102, -t81 * qJD(3) - t83 * t103, (t81 * t92 + t82 * t93) * qJD(1), t82 * qJD(3) - t83 * t104, t82 ^ 2 / 0.2e1 + t83 ^ 2 / 0.2e1 + t81 ^ 2 / 0.2e1, -t78 * qJD(3) + t79 * t103, t80 * qJD(3) + t79 * t104, (-t78 * t92 - t80 * t93) * qJD(1), t80 ^ 2 / 0.2e1 + t78 ^ 2 / 0.2e1 + t79 ^ 2 / 0.2e1;];
T_reg = t1;
