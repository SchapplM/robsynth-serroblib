% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x16]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:11
% EndTime: 2019-12-05 16:40:11
% DurationCPUTime: 0.10s
% Computational Cost: add. (54->18), mult. (92->46), div. (0->0), fcn. (37->4), ass. (0->18)
t72 = qJD(2) + qJD(3);
t71 = t72 ^ 2;
t84 = t71 / 0.2e1;
t81 = pkin(2) * qJD(2);
t77 = cos(qJ(3)) * t81;
t83 = (-t72 * pkin(3) - t77) * t72;
t78 = sin(qJ(3)) * t81;
t67 = t72 * pkin(7) + t78;
t73 = sin(qJ(4));
t75 = cos(qJ(4));
t82 = t73 * qJD(1) + t75 * t67;
t80 = qJ(5) * t72;
t79 = qJD(4) * t72;
t70 = t75 * qJD(1);
t65 = -t77 + qJD(5) + (-pkin(4) * t75 - pkin(3)) * t72;
t64 = t75 * t80 + t82;
t63 = qJD(4) * pkin(4) + t70 + (-t67 - t80) * t73;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, t84, t72 * t77, -t72 * t78, t73 ^ 2 * t84, t73 * t71 * t75, t73 * t79, t75 * t79, qJD(4) ^ 2 / 0.2e1, -t75 * t83 + (-t67 * t73 + t70) * qJD(4), -qJD(4) * t82 + t73 * t83, (-t63 * t73 + t64 * t75) * t72, t64 ^ 2 / 0.2e1 + t63 ^ 2 / 0.2e1 + t65 ^ 2 / 0.2e1;];
T_reg = t1;
