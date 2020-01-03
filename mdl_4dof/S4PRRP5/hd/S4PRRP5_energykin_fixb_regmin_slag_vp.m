% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% T_reg [1x13]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:21
% EndTime: 2019-12-31 16:29:21
% DurationCPUTime: 0.05s
% Computational Cost: add. (23->13), mult. (73->41), div. (0->0), fcn. (31->4), ass. (0->17)
t66 = qJD(2) ^ 2;
t73 = t66 / 0.2e1;
t65 = cos(qJ(2));
t70 = t65 * qJD(1);
t72 = qJD(2) * (-qJD(2) * pkin(2) - t70);
t63 = sin(qJ(2));
t60 = qJD(2) * pkin(5) + t63 * qJD(1);
t71 = qJD(3) * t60;
t69 = qJD(1) * qJD(2);
t68 = qJD(2) * qJD(3);
t67 = qJ(4) * qJD(2) + t60;
t64 = cos(qJ(3));
t62 = sin(qJ(3));
t59 = -t70 + qJD(4) + (-pkin(3) * t64 - pkin(2)) * qJD(2);
t58 = t67 * t64;
t57 = qJD(3) * pkin(3) - t62 * t67;
t1 = [qJD(1) ^ 2 / 0.2e1, t73, t65 * t69, -t63 * t69, t62 ^ 2 * t73, t62 * t66 * t64, t62 * t68, t64 * t68, qJD(3) ^ 2 / 0.2e1, -t62 * t71 - t64 * t72, t62 * t72 - t64 * t71, (-t57 * t62 + t58 * t64) * qJD(2), t58 ^ 2 / 0.2e1 + t57 ^ 2 / 0.2e1 + t59 ^ 2 / 0.2e1;];
T_reg = t1;
