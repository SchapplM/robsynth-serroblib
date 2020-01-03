% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPP4
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:17
% EndTime: 2019-12-31 16:59:17
% DurationCPUTime: 0.05s
% Computational Cost: add. (57->23), mult. (146->57), div. (0->0), fcn. (50->2), ass. (0->21)
t69 = qJD(1) ^ 2;
t80 = t69 / 0.2e1;
t79 = pkin(2) + pkin(3);
t68 = cos(qJ(2));
t78 = t68 * t69;
t76 = qJD(1) * t68;
t63 = pkin(5) * t76 + qJD(2) * qJ(3);
t67 = sin(qJ(2));
t77 = qJD(1) * t67;
t75 = pkin(5) * t77 + qJD(3);
t74 = qJ(4) * qJD(1);
t73 = qJD(1) * qJD(2);
t72 = t67 * t73;
t71 = t68 * t73;
t70 = qJ(3) * t67 + pkin(1);
t62 = -qJD(2) * pkin(2) + t75;
t61 = (-pkin(2) * t68 - t70) * qJD(1);
t60 = -t68 * t74 + t63;
t59 = -t79 * qJD(2) - t67 * t74 + t75;
t58 = qJD(4) + (t79 * t68 + t70) * qJD(1);
t1 = [t80, 0, 0, t67 ^ 2 * t80, t67 * t78, t72, t71, qJD(2) ^ 2 / 0.2e1, pkin(1) * t78 - pkin(5) * t72, -t69 * pkin(1) * t67 - pkin(5) * t71, -t62 * qJD(2) - t61 * t76, (t62 * t67 + t63 * t68) * qJD(1), t63 * qJD(2) - t61 * t77, t63 ^ 2 / 0.2e1 + t61 ^ 2 / 0.2e1 + t62 ^ 2 / 0.2e1, -t59 * qJD(2) + t58 * t76, t60 * qJD(2) + t58 * t77, (-t59 * t67 - t60 * t68) * qJD(1), t60 ^ 2 / 0.2e1 + t59 ^ 2 / 0.2e1 + t58 ^ 2 / 0.2e1;];
T_reg = t1;
