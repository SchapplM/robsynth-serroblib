% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:36
% EndTime: 2019-12-31 17:00:36
% DurationCPUTime: 0.07s
% Computational Cost: add. (57->24), mult. (146->56), div. (0->0), fcn. (50->2), ass. (0->21)
t66 = qJD(1) ^ 2;
t77 = t66 / 0.2e1;
t65 = cos(qJ(2));
t76 = t65 * t66;
t75 = -pkin(2) - qJ(4);
t64 = sin(qJ(2));
t74 = qJD(1) * t64;
t73 = qJD(1) * t65;
t72 = pkin(5) * t74 + qJD(3);
t71 = qJD(2) * qJ(3);
t70 = qJD(1) * qJD(2);
t69 = t64 * t70;
t68 = t65 * t70;
t67 = -qJ(3) * t64 - pkin(1);
t62 = -pkin(5) * t73 - t71;
t61 = -qJD(2) * pkin(2) + t72;
t60 = (-pkin(2) * t65 + t67) * qJD(1);
t59 = t71 + qJD(4) + (pkin(3) + pkin(5)) * t73;
t58 = pkin(3) * t74 + t75 * qJD(2) + t72;
t57 = (t75 * t65 + t67) * qJD(1);
t1 = [t77, 0, 0, t64 ^ 2 * t77, t64 * t76, t69, t68, qJD(2) ^ 2 / 0.2e1, pkin(1) * t76 - pkin(5) * t69, -t66 * pkin(1) * t64 - pkin(5) * t68, (t61 * t64 - t62 * t65) * qJD(1), t61 * qJD(2) + t60 * t73, -t62 * qJD(2) - t60 * t74, t60 ^ 2 / 0.2e1 + t62 ^ 2 / 0.2e1 + t61 ^ 2 / 0.2e1, (t58 * t64 + t59 * t65) * qJD(1), t59 * qJD(2) - t57 * t74, -t58 * qJD(2) - t57 * t73, t57 ^ 2 / 0.2e1 + t58 ^ 2 / 0.2e1 + t59 ^ 2 / 0.2e1;];
T_reg = t1;
