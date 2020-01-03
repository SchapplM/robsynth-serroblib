% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRRP6
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
% T_reg [1x15]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:42
% EndTime: 2019-12-31 16:30:42
% DurationCPUTime: 0.06s
% Computational Cost: add. (36->16), mult. (95->47), div. (0->0), fcn. (40->4), ass. (0->18)
t74 = qJD(2) ^ 2;
t81 = t74 / 0.2e1;
t70 = sin(qJ(3));
t80 = qJD(2) * t70;
t72 = cos(qJ(3));
t79 = qJD(2) * t72;
t71 = sin(qJ(2));
t68 = qJD(2) * pkin(5) + t71 * qJD(1);
t78 = qJD(3) * t68;
t73 = cos(qJ(2));
t77 = t73 * qJD(1);
t76 = qJD(1) * qJD(2);
t75 = qJD(2) * qJD(3);
t69 = -qJD(2) * pkin(2) - t77;
t67 = qJD(3) * qJ(4) + t72 * t68;
t66 = -qJD(3) * pkin(3) + t70 * t68 + qJD(4);
t65 = -t77 + (-pkin(3) * t72 - qJ(4) * t70 - pkin(2)) * qJD(2);
t1 = [qJD(1) ^ 2 / 0.2e1, t81, t73 * t76, -t71 * t76, t70 ^ 2 * t81, t70 * t74 * t72, t70 * t75, t72 * t75, qJD(3) ^ 2 / 0.2e1, -t69 * t79 - t70 * t78, t69 * t80 - t72 * t78, -t66 * qJD(3) - t65 * t79, (t66 * t70 + t67 * t72) * qJD(2), t67 * qJD(3) - t65 * t80, t67 ^ 2 / 0.2e1 + t65 ^ 2 / 0.2e1 + t66 ^ 2 / 0.2e1;];
T_reg = t1;
