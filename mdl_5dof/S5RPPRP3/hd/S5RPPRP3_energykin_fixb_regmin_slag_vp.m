% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:01
% EndTime: 2019-12-31 17:51:01
% DurationCPUTime: 0.06s
% Computational Cost: add. (53->24), mult. (122->53), div. (0->0), fcn. (41->4), ass. (0->21)
t80 = qJD(1) ^ 2;
t89 = t80 / 0.2e1;
t77 = cos(pkin(7));
t84 = -pkin(1) * t77 - pkin(2);
t70 = qJD(3) + (-pkin(6) + t84) * qJD(1);
t78 = sin(qJ(4));
t79 = cos(qJ(4));
t88 = t79 * qJD(2) + t78 * t70;
t76 = sin(pkin(7));
t83 = -pkin(1) * t76 - qJ(3);
t72 = t83 * qJD(1);
t87 = t72 * qJD(1);
t86 = qJ(5) * qJD(1);
t85 = qJD(1) * qJD(4);
t82 = -t78 * qJD(2) + t79 * t70;
t75 = qJD(2) ^ 2 / 0.2e1;
t71 = t84 * qJD(1) + qJD(3);
t69 = qJD(5) + (pkin(4) * t78 - t83) * qJD(1);
t66 = -t78 * t86 + t88;
t65 = qJD(4) * pkin(4) - t79 * t86 + t82;
t1 = [t89, 0, 0, t75 + (t76 ^ 2 / 0.2e1 + t77 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t80, t71 * qJD(1), -t87, t75 + t72 ^ 2 / 0.2e1 + t71 ^ 2 / 0.2e1, t79 ^ 2 * t89, -t79 * t80 * t78, t79 * t85, -t78 * t85, qJD(4) ^ 2 / 0.2e1, t82 * qJD(4) - t78 * t87, -t88 * qJD(4) - t79 * t87, (-t65 * t79 - t66 * t78) * qJD(1), t66 ^ 2 / 0.2e1 + t65 ^ 2 / 0.2e1 + t69 ^ 2 / 0.2e1;];
T_reg = t1;
