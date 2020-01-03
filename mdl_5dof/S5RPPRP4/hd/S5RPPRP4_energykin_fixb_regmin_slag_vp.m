% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:19
% EndTime: 2019-12-31 17:52:19
% DurationCPUTime: 0.06s
% Computational Cost: add. (76->27), mult. (148->58), div. (0->0), fcn. (52->4), ass. (0->22)
t88 = qJD(1) ^ 2;
t94 = t88 / 0.2e1;
t77 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t83 = sin(pkin(7));
t84 = cos(pkin(7));
t91 = qJ(2) * qJD(1);
t75 = t83 * t77 + t84 * t91;
t73 = -qJD(1) * pkin(6) + t75;
t86 = sin(qJ(4));
t87 = cos(qJ(4));
t93 = t86 * qJD(3) + t87 * t73;
t92 = qJD(1) * t87;
t90 = qJ(5) * qJD(1);
t89 = qJD(1) * qJD(4);
t74 = t84 * t77 - t83 * t91;
t72 = qJD(1) * pkin(3) - t74;
t82 = t87 * qJD(3);
t80 = -qJD(1) * pkin(1) + qJD(2);
t70 = pkin(4) * t92 + qJD(5) + t72;
t69 = -t87 * t90 + t93;
t68 = qJD(4) * pkin(4) + t82 + (-t73 + t90) * t86;
t1 = [t94, 0, 0, -t80 * qJD(1), t88 * qJ(2), qJ(2) ^ 2 * t94 + t80 ^ 2 / 0.2e1, -t74 * qJD(1), t75 * qJD(1), t75 ^ 2 / 0.2e1 + t74 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t86 ^ 2 * t94, t86 * t88 * t87, -t86 * t89, -t87 * t89, qJD(4) ^ 2 / 0.2e1, t72 * t92 + (-t86 * t73 + t82) * qJD(4), -t72 * t86 * qJD(1) - t93 * qJD(4), (t68 * t86 - t69 * t87) * qJD(1), t69 ^ 2 / 0.2e1 + t68 ^ 2 / 0.2e1 + t70 ^ 2 / 0.2e1;];
T_reg = t1;
