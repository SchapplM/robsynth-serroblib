% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% T_reg [1x12]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PPRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:39
% EndTime: 2019-12-31 16:18:39
% DurationCPUTime: 0.04s
% Computational Cost: add. (17->12), mult. (67->39), div. (0->0), fcn. (37->6), ass. (0->15)
t55 = sin(pkin(7));
t56 = cos(pkin(7));
t58 = sin(qJ(3));
t60 = cos(qJ(3));
t68 = (t55 * t58 - t56 * t60) * qJD(1);
t61 = qJD(3) ^ 2;
t67 = t61 / 0.2e1;
t66 = (t55 * t60 + t56 * t58) * qJD(1);
t65 = (-qJD(3) * pkin(3) + t68) * qJD(3);
t64 = qJD(3) * qJD(4);
t62 = qJD(1) ^ 2;
t59 = cos(qJ(4));
t57 = sin(qJ(4));
t52 = qJD(3) * pkin(5) + t66;
t1 = [t62 / 0.2e1, qJD(2) ^ 2 / 0.2e1 + (t55 ^ 2 / 0.2e1 + t56 ^ 2 / 0.2e1) * t62, t67, -qJD(3) * t68, -t66 * qJD(3), t57 ^ 2 * t67, t57 * t61 * t59, t57 * t64, t59 * t64, qJD(4) ^ 2 / 0.2e1, (t59 * qJD(2) - t57 * t52) * qJD(4) - t59 * t65, -(t57 * qJD(2) + t59 * t52) * qJD(4) + t57 * t65;];
T_reg = t1;
