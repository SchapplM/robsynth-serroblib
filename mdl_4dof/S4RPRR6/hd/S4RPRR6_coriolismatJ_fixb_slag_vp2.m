% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR6
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
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:30
% EndTime: 2019-12-31 16:52:31
% DurationCPUTime: 0.56s
% Computational Cost: add. (2136->82), mult. (4224->118), div. (0->0), fcn. (4704->6), ass. (0->53)
t74 = sin(pkin(7));
t95 = pkin(5) + qJ(2);
t70 = t95 * t74;
t75 = cos(pkin(7));
t71 = t95 * t75;
t77 = sin(qJ(3));
t79 = cos(qJ(3));
t61 = -t77 * t70 + t71 * t79;
t66 = -t77 * t74 + t75 * t79;
t43 = t66 * pkin(6) + t61;
t76 = sin(qJ(4));
t78 = cos(qJ(4));
t60 = -t79 * t70 - t71 * t77;
t67 = t74 * t79 + t77 * t75;
t80 = -pkin(6) * t67 + t60;
t118 = -t43 * t76 + t78 * t80;
t34 = t43 * t78 + t76 * t80;
t59 = t66 * t76 + t67 * t78;
t83 = t78 * t66 - t67 * t76;
t3 = -t34 * mrSges(5,1) - t118 * mrSges(5,2) + Ifges(5,5) * t83 - Ifges(5,6) * t59;
t128 = t3 * qJD(4);
t86 = -pkin(2) * t75 - pkin(1);
t81 = pkin(3) * t66 - t86;
t119 = t83 * mrSges(5,2);
t84 = t59 * mrSges(5,1) + t119;
t126 = t81 * t84;
t124 = (-Ifges(5,1) + Ifges(5,2)) * t83;
t115 = t59 / 0.2e1;
t87 = t119 / 0.2e1;
t99 = t83 ^ 2;
t117 = t66 ^ 2;
t116 = m(5) * pkin(3);
t102 = t59 ^ 2;
t101 = t59 * t78;
t98 = t83 * t76;
t5 = (t99 + t102) * mrSges(5,3) + (t67 ^ 2 + t117) * mrSges(4,3) + m(5) * (-t118 * t59 + t34 * t83) + m(4) * (-t60 * t67 + t61 * t66) + (m(3) * qJ(2) + mrSges(3,3)) * (t74 ^ 2 + t75 ^ 2);
t93 = qJD(1) * t5;
t85 = t67 * mrSges(4,1) + t66 * mrSges(4,2);
t1 = t126 - Ifges(4,4) * t117 - t86 * t85 + (-t99 + t102) * Ifges(5,4) + (pkin(3) * mrSges(5,1) * t83 + t81 * t116 + (-Ifges(4,1) + Ifges(4,2)) * t66 + Ifges(4,4) * t67) * t67 + (-pkin(3) * mrSges(5,2) * t67 + t124) * t59;
t92 = t1 * qJD(1);
t2 = -Ifges(5,4) * t99 + t126 + (Ifges(5,4) * t59 + t124) * t59;
t91 = t2 * qJD(1);
t8 = (-t101 / 0.2e1 + t98 / 0.2e1 - t67 / 0.2e1) * t116 - t84 - t85;
t90 = t8 * qJD(1);
t10 = 0.2e1 * t115 * mrSges(5,1) + 0.2e1 * t87;
t89 = t10 * qJD(1);
t4 = (-t59 / 0.2e1 + t115) * Ifges(5,6);
t69 = (mrSges(5,1) * t76 + mrSges(5,2) * t78) * pkin(3);
t82 = -qJD(1) * t4 + qJD(3) * t69;
t68 = t69 * qJD(4);
t24 = (t67 + t98 - t101) * t116 / 0.2e1;
t11 = t87 - t119 / 0.2e1;
t6 = [qJD(2) * t5 - qJD(3) * t1 - qJD(4) * t2, qJD(3) * t24 + qJD(4) * t11 + t93, t24 * qJD(2) + t128 - t92 + (-t61 * mrSges(4,1) - t60 * mrSges(4,2) + Ifges(4,5) * t66 - Ifges(4,6) * t67 + (m(5) * (t118 * t76 - t34 * t78) + (-t59 * t76 - t78 * t83) * mrSges(5,3)) * pkin(3) + t3) * qJD(3), t11 * qJD(2) + t3 * qJD(3) + t128 - t91; -qJD(3) * t8 + qJD(4) * t10 - t93, 0, -t90, t89; qJD(2) * t8 + qJD(4) * t4 + t92, t90, -t68, -t68 - t82; -qJD(2) * t10 - qJD(3) * t4 + t91, -t89, t82, 0;];
Cq = t6;
