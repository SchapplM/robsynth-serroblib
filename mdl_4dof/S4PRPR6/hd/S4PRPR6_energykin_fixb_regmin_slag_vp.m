% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:40
% EndTime: 2019-12-31 16:24:41
% DurationCPUTime: 0.04s
% Computational Cost: add. (40->19), mult. (113->50), div. (0->0), fcn. (63->6), ass. (0->21)
t72 = sin(pkin(7));
t82 = qJD(2) * t72;
t73 = cos(pkin(7));
t81 = qJD(2) * t73;
t80 = qJD(1) * qJD(2);
t75 = sin(qJ(2));
t68 = qJD(2) * qJ(3) + t75 * qJD(1);
t79 = pkin(5) * qJD(2) + t68;
t77 = cos(qJ(2));
t78 = -t77 * qJD(1) + qJD(3);
t76 = cos(qJ(4));
t74 = sin(qJ(4));
t71 = t73 ^ 2;
t70 = t72 ^ 2;
t66 = -qJD(2) * pkin(2) + t78;
t65 = (-pkin(3) * t73 - pkin(2)) * qJD(2) + t78;
t64 = (t72 * t76 + t73 * t74) * qJD(2);
t63 = t74 * t82 - t76 * t81;
t62 = t79 * t73;
t61 = t79 * t72;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t77 * t80, -t75 * t80, -t66 * t81, t66 * t82, (t70 + t71) * t68 * qJD(2), t66 ^ 2 / 0.2e1 + (t71 / 0.2e1 + t70 / 0.2e1) * t68 ^ 2, t64 ^ 2 / 0.2e1, -t64 * t63, t64 * qJD(4), -t63 * qJD(4), qJD(4) ^ 2 / 0.2e1, (-t76 * t61 - t74 * t62) * qJD(4) + t65 * t63, -(-t74 * t61 + t76 * t62) * qJD(4) + t65 * t64;];
T_reg = t1;
