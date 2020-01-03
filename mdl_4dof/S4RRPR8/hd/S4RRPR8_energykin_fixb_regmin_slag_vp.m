% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:22
% EndTime: 2019-12-31 17:08:22
% DurationCPUTime: 0.05s
% Computational Cost: add. (56->26), mult. (155->64), div. (0->0), fcn. (71->4), ass. (0->25)
t93 = qJD(1) ^ 2;
t103 = t93 / 0.2e1;
t102 = pkin(2) + pkin(3);
t92 = cos(qJ(2));
t101 = t92 * t93;
t99 = qJD(1) * t92;
t83 = pkin(5) * t99 + qJD(2) * qJ(3);
t90 = sin(qJ(2));
t100 = qJD(1) * t90;
t98 = pkin(5) * t100 + qJD(3);
t97 = qJD(1) * qJD(2);
t96 = t90 * t97;
t95 = t92 * t97;
t94 = qJ(3) * t90 + pkin(1);
t91 = cos(qJ(4));
t89 = sin(qJ(4));
t86 = qJD(2) - qJD(4);
t82 = -qJD(2) * pkin(2) + t98;
t81 = (-pkin(2) * t92 - t94) * qJD(1);
t80 = -pkin(6) * t99 + t83;
t79 = (-t89 * t92 + t90 * t91) * qJD(1);
t78 = (t89 * t90 + t91 * t92) * qJD(1);
t77 = -pkin(6) * t100 - t102 * qJD(2) + t98;
t76 = (t102 * t92 + t94) * qJD(1);
t1 = [t103, 0, 0, t90 ^ 2 * t103, t90 * t101, t96, t95, qJD(2) ^ 2 / 0.2e1, pkin(1) * t101 - pkin(5) * t96, -t93 * pkin(1) * t90 - pkin(5) * t95, -t82 * qJD(2) - t81 * t99, (t82 * t90 + t83 * t92) * qJD(1), t83 * qJD(2) - t81 * t100, t83 ^ 2 / 0.2e1 + t81 ^ 2 / 0.2e1 + t82 ^ 2 / 0.2e1, t79 ^ 2 / 0.2e1, -t79 * t78, -t79 * t86, t78 * t86, t86 ^ 2 / 0.2e1, t76 * t78 - (t91 * t77 - t89 * t80) * t86, t76 * t79 + (t89 * t77 + t91 * t80) * t86;];
T_reg = t1;
