% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% T_reg [1x18]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:01
% EndTime: 2019-12-31 16:45:01
% DurationCPUTime: 0.07s
% Computational Cost: add. (80->26), mult. (229->58), div. (0->0), fcn. (127->4), ass. (0->22)
t87 = pkin(5) + qJ(2);
t77 = sin(pkin(6));
t85 = qJD(1) * t77;
t70 = t87 * t85;
t78 = cos(pkin(6));
t84 = qJD(1) * t78;
t71 = t87 * t84;
t79 = sin(qJ(3));
t80 = cos(qJ(3));
t86 = -t79 * t70 + t80 * t71;
t83 = -t80 * t70 - t79 * t71;
t72 = qJD(2) + (-pkin(2) * t78 - pkin(1)) * qJD(1);
t81 = qJD(1) ^ 2;
t76 = t78 ^ 2;
t75 = t77 ^ 2;
t74 = -qJD(1) * pkin(1) + qJD(2);
t69 = (t77 * t80 + t78 * t79) * qJD(1);
t68 = t79 * t85 - t80 * t84;
t65 = qJD(3) * qJ(4) + t86;
t64 = -qJD(3) * pkin(3) + qJD(4) - t83;
t63 = t68 * pkin(3) - t69 * qJ(4) + t72;
t1 = [t81 / 0.2e1, 0, 0, -t74 * t84, t74 * t85, (t75 + t76) * t81 * qJ(2), t74 ^ 2 / 0.2e1 + (t76 / 0.2e1 + t75 / 0.2e1) * qJ(2) ^ 2 * t81, t69 ^ 2 / 0.2e1, -t69 * t68, t69 * qJD(3), -t68 * qJD(3), qJD(3) ^ 2 / 0.2e1, qJD(3) * t83 + t72 * t68, -t86 * qJD(3) + t72 * t69, -t64 * qJD(3) + t63 * t68, t64 * t69 - t65 * t68, t65 * qJD(3) - t63 * t69, t65 ^ 2 / 0.2e1 + t63 ^ 2 / 0.2e1 + t64 ^ 2 / 0.2e1;];
T_reg = t1;
