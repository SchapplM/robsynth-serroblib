% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRPRR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:44
% EndTime: 2019-12-31 17:39:44
% DurationCPUTime: 0.08s
% Computational Cost: add. (45->16), mult. (78->40), div. (0->0), fcn. (23->4), ass. (0->19)
t66 = -qJD(2) + qJD(4);
t65 = t66 ^ 2;
t79 = t65 / 0.2e1;
t72 = qJD(2) ^ 2;
t78 = t72 / 0.2e1;
t62 = qJD(3) + (-pkin(2) - pkin(3)) * qJD(2);
t69 = sin(qJ(4));
t71 = cos(qJ(4));
t74 = qJ(3) * qJD(2);
t73 = t71 * t62 - t69 * t74;
t77 = (-t66 * pkin(4) - t73) * t66;
t76 = t69 * t62 + t71 * t74;
t75 = qJD(5) * t66;
t70 = cos(qJ(5));
t68 = sin(qJ(5));
t67 = qJD(1) ^ 2 / 0.2e1;
t64 = -qJD(2) * pkin(2) + qJD(3);
t60 = t66 * pkin(7) + t76;
t1 = [t67, t78, 0, 0, -t64 * qJD(2), t72 * qJ(3), qJ(3) ^ 2 * t78 + t67 + t64 ^ 2 / 0.2e1, t79, t73 * t66, -t76 * t66, t68 ^ 2 * t79, t68 * t65 * t70, t68 * t75, t70 * t75, qJD(5) ^ 2 / 0.2e1, -t70 * t77 + (-qJD(1) * t70 - t60 * t68) * qJD(5), t68 * t77 - (-qJD(1) * t68 + t60 * t70) * qJD(5);];
T_reg = t1;
