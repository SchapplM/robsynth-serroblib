% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% T_reg [1x14]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PPRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:40
% EndTime: 2019-12-31 17:34:41
% DurationCPUTime: 0.06s
% Computational Cost: add. (30->18), mult. (83->45), div. (0->0), fcn. (37->4), ass. (0->19)
t73 = qJD(3) ^ 2;
t81 = t73 / 0.2e1;
t71 = cos(qJ(4));
t80 = qJD(1) * t71;
t72 = cos(qJ(3));
t79 = qJD(2) * t72;
t78 = qJD(3) * (-qJD(3) * pkin(3) - t79);
t77 = qJ(5) * qJD(3);
t76 = qJD(2) * qJD(3);
t75 = qJD(3) * qJD(4);
t70 = sin(qJ(3));
t66 = qJD(3) * pkin(6) + qJD(2) * t70;
t69 = sin(qJ(4));
t74 = -qJD(1) * t69 + t71 * t66;
t68 = qJD(1) ^ 2 / 0.2e1;
t64 = -t79 + qJD(5) + (-pkin(4) * t71 - pkin(3)) * qJD(3);
t63 = t71 * t77 + t74;
t62 = qJD(4) * pkin(4) - t80 + (-t66 - t77) * t69;
t1 = [t68, t68 + qJD(2) ^ 2 / 0.2e1, t81, t72 * t76, -t70 * t76, t69 ^ 2 * t81, t69 * t73 * t71, t69 * t75, t71 * t75, qJD(4) ^ 2 / 0.2e1, -t71 * t78 + (-t66 * t69 - t80) * qJD(4), -t74 * qJD(4) + t69 * t78, (-t62 * t69 + t63 * t71) * qJD(3), t63 ^ 2 / 0.2e1 + t62 ^ 2 / 0.2e1 + t64 ^ 2 / 0.2e1;];
T_reg = t1;
