% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRPR5
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
% T_reg [1x12]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:17
% EndTime: 2019-12-31 16:23:17
% DurationCPUTime: 0.05s
% Computational Cost: add. (23->12), mult. (71->40), div. (0->0), fcn. (37->6), ass. (0->17)
t70 = qJD(2) ^ 2;
t75 = t70 / 0.2e1;
t69 = cos(qJ(2));
t62 = qJD(2) * pkin(2) + t69 * qJD(1);
t64 = sin(pkin(7));
t65 = cos(pkin(7));
t67 = sin(qJ(2));
t74 = qJD(1) * t67;
t60 = t64 * t62 + t65 * t74;
t59 = t65 * t62 - t64 * t74;
t73 = qJD(2) * (-qJD(2) * pkin(3) - t59);
t72 = qJD(1) * qJD(2);
t71 = qJD(2) * qJD(4);
t68 = cos(qJ(4));
t66 = sin(qJ(4));
t58 = qJD(2) * pkin(5) + t60;
t1 = [qJD(1) ^ 2 / 0.2e1, t75, t69 * t72, -t67 * t72, t60 ^ 2 / 0.2e1 + t59 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t66 ^ 2 * t75, t66 * t70 * t68, t66 * t71, t68 * t71, qJD(4) ^ 2 / 0.2e1, (t68 * qJD(3) - t66 * t58) * qJD(4) - t68 * t73, -(t66 * qJD(3) + t68 * t58) * qJD(4) + t66 * t73;];
T_reg = t1;
