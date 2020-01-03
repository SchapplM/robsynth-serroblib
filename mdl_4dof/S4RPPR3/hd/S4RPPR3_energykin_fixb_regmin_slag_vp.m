% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% 
% Output:
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:56
% EndTime: 2019-12-31 16:37:56
% DurationCPUTime: 0.04s
% Computational Cost: add. (49->24), mult. (140->57), div. (0->0), fcn. (73->6), ass. (0->21)
t72 = sin(pkin(6));
t67 = (pkin(1) * t72 + qJ(3)) * qJD(1);
t71 = sin(pkin(7));
t73 = cos(pkin(7));
t61 = t71 * qJD(2) + t73 * t67;
t81 = qJD(1) * t71;
t80 = qJD(1) * t73;
t74 = cos(pkin(6));
t79 = -pkin(1) * t74 - pkin(2);
t77 = qJD(1) ^ 2;
t76 = cos(qJ(4));
t75 = sin(qJ(4));
t70 = t73 * qJD(2);
t66 = t79 * qJD(1) + qJD(3);
t64 = (t71 * t76 + t73 * t75) * qJD(1);
t63 = t75 * t81 - t76 * t80;
t62 = qJD(3) + (-pkin(3) * t73 + t79) * qJD(1);
t60 = -t71 * t67 + t70;
t59 = pkin(5) * t80 + t61;
t58 = t70 + (-pkin(5) * qJD(1) - t67) * t71;
t1 = [t77 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t72 ^ 2 / 0.2e1 + t74 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t77, -t66 * t80, t66 * t81, (-t60 * t71 + t61 * t73) * qJD(1), t61 ^ 2 / 0.2e1 + t60 ^ 2 / 0.2e1 + t66 ^ 2 / 0.2e1, t64 ^ 2 / 0.2e1, -t64 * t63, t64 * qJD(4), -t63 * qJD(4), qJD(4) ^ 2 / 0.2e1, t62 * t63 + (t76 * t58 - t75 * t59) * qJD(4), t62 * t64 - (t75 * t58 + t76 * t59) * qJD(4);];
T_reg = t1;
