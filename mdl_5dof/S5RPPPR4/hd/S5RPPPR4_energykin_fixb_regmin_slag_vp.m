% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:15
% EndTime: 2019-12-31 17:45:15
% DurationCPUTime: 0.06s
% Computational Cost: add. (73->27), mult. (168->62), div. (0->0), fcn. (77->6), ass. (0->23)
t86 = sin(pkin(7));
t80 = (-pkin(1) * t86 - qJ(3)) * qJD(1);
t88 = cos(pkin(7));
t93 = -pkin(1) * t88 - pkin(2);
t77 = qJD(3) + (-qJ(4) + t93) * qJD(1);
t85 = sin(pkin(8));
t87 = cos(pkin(8));
t71 = t87 * qJD(2) + t85 * t77;
t95 = qJD(1) * t85;
t94 = qJD(1) * t87;
t78 = qJD(4) - t80;
t70 = -t85 * qJD(2) + t87 * t77;
t91 = qJD(1) ^ 2;
t90 = cos(qJ(5));
t89 = sin(qJ(5));
t83 = qJD(2) ^ 2 / 0.2e1;
t79 = t93 * qJD(1) + qJD(3);
t76 = (-t85 * t89 + t87 * t90) * qJD(1);
t75 = (t85 * t90 + t87 * t89) * qJD(1);
t74 = pkin(4) * t95 + t78;
t69 = -pkin(6) * t95 + t71;
t68 = -pkin(6) * t94 + t70;
t1 = [t91 / 0.2e1, 0, 0, t83 + (t86 ^ 2 / 0.2e1 + t88 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t91, t79 * qJD(1), -t80 * qJD(1), t83 + t80 ^ 2 / 0.2e1 + t79 ^ 2 / 0.2e1, t78 * t95, t78 * t94, (-t70 * t87 - t71 * t85) * qJD(1), t71 ^ 2 / 0.2e1 + t70 ^ 2 / 0.2e1 + t78 ^ 2 / 0.2e1, t76 ^ 2 / 0.2e1, -t76 * t75, t76 * qJD(5), -t75 * qJD(5), qJD(5) ^ 2 / 0.2e1, t74 * t75 + (t90 * t68 - t89 * t69) * qJD(5), t74 * t76 - (t89 * t68 + t90 * t69) * qJD(5);];
T_reg = t1;
