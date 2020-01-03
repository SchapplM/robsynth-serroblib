% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x14]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:12
% EndTime: 2019-12-31 16:48:12
% DurationCPUTime: 0.06s
% Computational Cost: add. (36->13), mult. (89->42), div. (0->0), fcn. (37->6), ass. (0->18)
t60 = qJD(1) + qJD(3);
t59 = t60 ^ 2;
t74 = t59 / 0.2e1;
t62 = cos(pkin(7));
t57 = (pkin(1) * t62 + pkin(2)) * qJD(1);
t64 = sin(qJ(3));
t66 = cos(qJ(3));
t61 = sin(pkin(7));
t70 = pkin(1) * qJD(1) * t61;
t69 = t66 * t57 - t64 * t70;
t73 = (-t60 * pkin(3) - t69) * t60;
t72 = t64 * t57 + t66 * t70;
t71 = qJD(4) * t60;
t67 = qJD(1) ^ 2;
t65 = cos(qJ(4));
t63 = sin(qJ(4));
t55 = t60 * pkin(6) + t72;
t1 = [t67 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t61 ^ 2 / 0.2e1 + t62 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t67, t74, t69 * t60, -t72 * t60, t63 ^ 2 * t74, t63 * t59 * t65, t63 * t71, t65 * t71, qJD(4) ^ 2 / 0.2e1, -t65 * t73 + (t65 * qJD(2) - t63 * t55) * qJD(4), t63 * t73 - (t63 * qJD(2) + t65 * t55) * qJD(4);];
T_reg = t1;
