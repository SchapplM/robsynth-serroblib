% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:11
% EndTime: 2019-12-31 17:15:13
% DurationCPUTime: 0.70s
% Computational Cost: add. (1735->104), mult. (3522->138), div. (0->0), fcn. (3332->4), ass. (0->58)
t111 = sin(qJ(2));
t112 = cos(qJ(2));
t74 = sin(qJ(3));
t75 = cos(qJ(3));
t62 = -t74 * t111 + t75 * t112;
t72 = -t112 * pkin(2) - pkin(1);
t46 = -t62 * pkin(3) + t72;
t63 = -t75 * t111 - t74 * t112;
t131 = m(5) * t46 - mrSges(5,1) * t62 - mrSges(5,2) * t63;
t86 = t111 * pkin(5);
t67 = -t111 * pkin(6) - t86;
t88 = t112 * pkin(5);
t68 = t112 * pkin(6) + t88;
t117 = t75 * t67 - t74 * t68;
t120 = t63 * qJ(4) + t117;
t42 = t74 * t67 + t75 * t68;
t31 = t62 * qJ(4) + t42;
t130 = (Ifges(4,6) + Ifges(5,6)) * t63 + (Ifges(4,5) + Ifges(5,5)) * t62 - t117 * mrSges(4,2) - t120 * mrSges(5,2) - t42 * mrSges(4,1) - t31 * mrSges(5,1);
t110 = mrSges(5,3) * t62;
t124 = -m(5) * t31 - t110;
t122 = Ifges(5,4) + Ifges(4,4);
t123 = mrSges(5,3) * t120 + t122 * t62;
t113 = t75 * pkin(2);
t119 = pkin(2) * t74;
t116 = m(5) / 0.2e1;
t114 = pkin(3) * t31;
t101 = t63 * t74;
t71 = pkin(3) + t113;
t100 = t71 * t31;
t99 = t71 * t63;
t97 = mrSges(4,2) + mrSges(5,2);
t6 = m(5) * (t120 * t63 + t31 * t62) + (t62 ^ 2 + t63 ^ 2) * mrSges(5,3);
t94 = qJD(1) * t6;
t87 = t111 * pkin(2);
t52 = -t63 * pkin(3) + t87;
t55 = t62 * mrSges(5,2);
t32 = -t63 * mrSges(5,1) + t55;
t77 = -t120 * t110 + t46 * t32 + t72 * (-mrSges(4,1) * t63 + mrSges(4,2) * t62);
t83 = -Ifges(4,1) - Ifges(5,1) + Ifges(4,2) + Ifges(5,2);
t84 = t122 * t63;
t1 = -pkin(1) * (t111 * mrSges(3,1) + t112 * mrSges(3,2)) + m(4) * t72 * t87 + (-mrSges(4,1) * t87 + t123) * t62 + (-mrSges(4,2) * t87 + t83 * t62 - t84) * t63 + t77 + (-Ifges(3,2) + Ifges(3,1)) * t112 * t111 + (-t111 ^ 2 + t112 ^ 2) * Ifges(3,4) + t131 * t52;
t93 = t1 * qJD(1);
t2 = (-t131 * pkin(3) - t84) * t63 + (t83 * t63 + t123) * t62 + t77;
t92 = t2 * qJD(1);
t89 = t62 * t119;
t7 = 0.2e1 * (t89 / 0.4e1 + t99 / 0.4e1 - t52 / 0.4e1) * m(5) - t32;
t91 = t7 * qJD(1);
t85 = m(5) * pkin(3) + mrSges(5,1);
t11 = t85 * t63 - t55;
t90 = t11 * qJD(1);
t81 = -t71 / 0.2e1 + t113 / 0.2e1;
t21 = t97 * t113 + (mrSges(4,1) + mrSges(5,1) + (t71 - t113) * m(5)) * t119;
t76 = t31 * t113;
t4 = (-t100 / 0.2e1 + t76 / 0.2e1 + t114 / 0.2e1) * m(5) + (pkin(3) / 0.2e1 + t81) * t110;
t78 = -t4 * qJD(1) + t21 * qJD(2);
t9 = (t89 + t99 + t52) * t116;
t3 = (t76 - t100 - t114) * t116 + (-pkin(3) / 0.2e1 + t81) * t110 + t130;
t5 = [qJD(2) * t1 + qJD(3) * t2 + qJD(4) * t6, t3 * qJD(3) + t9 * qJD(4) + t93 + (-mrSges(3,1) * t88 + mrSges(3,2) * t86 + Ifges(3,5) * t112 - Ifges(3,6) * t111 + (mrSges(5,3) * t101 + (-t62 * t75 + t101) * mrSges(4,3) + m(5) * t120 * t74 + m(4) * (t117 * t74 - t42 * t75)) * pkin(2) + t124 * t71 + t130) * qJD(2), t92 + t3 * qJD(2) + (t124 * pkin(3) + t130) * qJD(3), qJD(2) * t9 + t94; qJD(3) * t4 + qJD(4) * t7 - t93, -t21 * qJD(3), (-t97 * t75 + (-mrSges(4,1) - t85) * t74) * qJD(3) * pkin(2) - t78, t91; -qJD(2) * t4 + qJD(4) * t11 - t92, t78, 0, t90; -qJD(2) * t7 - qJD(3) * t11 - t94, -t91, -t90, 0;];
Cq = t5;
