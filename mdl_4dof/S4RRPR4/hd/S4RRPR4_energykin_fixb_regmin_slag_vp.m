% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:33
% EndTime: 2019-12-31 17:02:33
% DurationCPUTime: 0.06s
% Computational Cost: add. (80->20), mult. (128->50), div. (0->0), fcn. (63->6), ass. (0->22)
t70 = qJD(1) + qJD(2);
t71 = sin(pkin(7));
t83 = t70 * t71;
t72 = cos(pkin(7));
t82 = t70 * t72;
t81 = pkin(1) * qJD(1);
t80 = sin(qJ(2)) * t81;
t79 = cos(qJ(2)) * t81;
t66 = qJ(3) * t70 + t80;
t78 = pkin(6) * t70 + t66;
t77 = qJD(3) - t79;
t75 = cos(qJ(4));
t73 = sin(qJ(4));
t69 = t72 ^ 2;
t68 = t71 ^ 2;
t64 = -pkin(2) * t70 + t77;
t63 = (t71 * t75 + t72 * t73) * t70;
t62 = t73 * t83 - t75 * t82;
t61 = (-pkin(3) * t72 - pkin(2)) * t70 + t77;
t60 = t78 * t72;
t59 = t78 * t71;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t70 ^ 2 / 0.2e1, t70 * t79, -t70 * t80, -t64 * t82, t64 * t83, (t68 + t69) * t70 * t66, t64 ^ 2 / 0.2e1 + (t69 / 0.2e1 + t68 / 0.2e1) * t66 ^ 2, t63 ^ 2 / 0.2e1, -t63 * t62, t63 * qJD(4), -t62 * qJD(4), qJD(4) ^ 2 / 0.2e1, t61 * t62 + (-t59 * t75 - t60 * t73) * qJD(4), t61 * t63 - (-t59 * t73 + t60 * t75) * qJD(4);];
T_reg = t1;
