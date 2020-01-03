% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:13
% EndTime: 2019-12-31 17:22:13
% DurationCPUTime: 0.07s
% Computational Cost: add. (55->11), mult. (81->35), div. (0->0), fcn. (35->6), ass. (0->19)
t61 = qJD(1) + qJD(2);
t60 = qJD(3) + t61;
t59 = t60 ^ 2;
t76 = t59 / 0.2e1;
t73 = pkin(1) * qJD(1);
t69 = cos(qJ(2)) * t73;
t57 = pkin(2) * t61 + t69;
t63 = sin(qJ(3));
t66 = cos(qJ(3));
t70 = sin(qJ(2)) * t73;
t68 = t57 * t66 - t63 * t70;
t75 = (-pkin(3) * t60 - t68) * t60;
t74 = t63 * t57 + t66 * t70;
t62 = sin(qJ(4));
t72 = qJD(4) * t62;
t65 = cos(qJ(4));
t71 = qJD(4) * t65;
t55 = pkin(7) * t60 + t74;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t61 ^ 2 / 0.2e1, t61 * t69, -t61 * t70, t76, t68 * t60, -t74 * t60, t62 ^ 2 * t76, t62 * t59 * t65, t60 * t72, t60 * t71, qJD(4) ^ 2 / 0.2e1, -t55 * t72 - t65 * t75, -t55 * t71 + t62 * t75;];
T_reg = t1;
