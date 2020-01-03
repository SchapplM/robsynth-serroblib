% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:11
% EndTime: 2019-12-31 17:41:11
% DurationCPUTime: 0.05s
% Computational Cost: add. (71->25), mult. (162->61), div. (0->0), fcn. (64->2), ass. (0->21)
t79 = qJD(2) ^ 2;
t90 = t79 / 0.2e1;
t89 = pkin(3) + pkin(4);
t78 = cos(qJ(3));
t88 = t78 * t79;
t77 = sin(qJ(3));
t85 = qJD(2) * t78;
t87 = pkin(6) * t85 + t77 * qJD(1);
t86 = qJD(2) * t77;
t84 = qJ(5) * qJD(2);
t83 = qJD(2) * qJD(3);
t72 = qJD(3) * qJ(4) + t87;
t82 = qJ(4) * t77 + pkin(2);
t81 = -pkin(6) * t86 + t78 * qJD(1);
t80 = qJD(4) - t81;
t71 = (-pkin(3) * t78 - t82) * qJD(2);
t70 = -qJD(3) * pkin(3) + t80;
t69 = -t78 * t84 + t72;
t68 = qJD(5) + (t89 * t78 + t82) * qJD(2);
t67 = -t89 * qJD(3) - t77 * t84 + t80;
t1 = [qJD(1) ^ 2 / 0.2e1, t90, 0, 0, t77 ^ 2 * t90, t77 * t88, t77 * t83, t78 * t83, qJD(3) ^ 2 / 0.2e1, pkin(2) * t88 + t81 * qJD(3), -t79 * pkin(2) * t77 - t87 * qJD(3), -t70 * qJD(3) - t71 * t85, (t70 * t77 + t72 * t78) * qJD(2), t72 * qJD(3) - t71 * t86, t72 ^ 2 / 0.2e1 + t71 ^ 2 / 0.2e1 + t70 ^ 2 / 0.2e1, -t67 * qJD(3) + t68 * t85, t69 * qJD(3) + t68 * t86, (-t67 * t77 - t69 * t78) * qJD(2), t69 ^ 2 / 0.2e1 + t67 ^ 2 / 0.2e1 + t68 ^ 2 / 0.2e1;];
T_reg = t1;
