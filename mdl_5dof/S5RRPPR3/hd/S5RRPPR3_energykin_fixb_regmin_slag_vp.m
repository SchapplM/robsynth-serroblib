% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:35
% EndTime: 2019-12-31 19:26:35
% DurationCPUTime: 0.09s
% Computational Cost: add. (77->20), mult. (127->45), div. (0->0), fcn. (53->6), ass. (0->21)
t70 = qJD(1) + qJD(2);
t69 = t70 ^ 2;
t84 = t69 / 0.2e1;
t82 = pkin(1) * qJD(1);
t79 = cos(qJ(2)) * t82;
t67 = t70 * pkin(2) + t79;
t72 = sin(pkin(8));
t73 = cos(pkin(8));
t80 = sin(qJ(2)) * t82;
t66 = t72 * t67 + t73 * t80;
t63 = t70 * qJ(4) + t66;
t83 = t63 * t70;
t81 = qJD(5) * t70;
t65 = t73 * t67 - t72 * t80;
t78 = qJD(4) - t65;
t76 = cos(qJ(5));
t74 = sin(qJ(5));
t71 = qJD(3) ^ 2 / 0.2e1;
t62 = -t70 * pkin(3) + t78;
t61 = (-pkin(3) - pkin(7)) * t70 + t78;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t84, t70 * t79, -t70 * t80, t66 ^ 2 / 0.2e1 + t65 ^ 2 / 0.2e1 + t71, t62 * t70, t83, t71 + t63 ^ 2 / 0.2e1 + t62 ^ 2 / 0.2e1, t76 ^ 2 * t84, -t76 * t69 * t74, t76 * t81, -t74 * t81, qJD(5) ^ 2 / 0.2e1, t74 * t83 + (-t74 * qJD(3) + t76 * t61) * qJD(5), t76 * t83 - (t76 * qJD(3) + t74 * t61) * qJD(5);];
T_reg = t1;
