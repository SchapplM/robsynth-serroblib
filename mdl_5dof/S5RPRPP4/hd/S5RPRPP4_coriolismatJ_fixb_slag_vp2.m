% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:19
% EndTime: 2019-12-31 18:14:20
% DurationCPUTime: 0.64s
% Computational Cost: add. (1638->102), mult. (2944->127), div. (0->0), fcn. (2794->4), ass. (0->59)
t103 = cos(qJ(3));
t66 = sin(pkin(7));
t67 = sin(qJ(3));
t93 = cos(pkin(7));
t55 = -t93 * t103 + t66 * t67;
t56 = t66 * t103 + t93 * t67;
t64 = t67 * pkin(3) + qJ(2);
t25 = pkin(4) * t56 + qJ(5) * t55 + t64;
t106 = m(6) * t25;
t31 = t56 * mrSges(6,1) + t55 * mrSges(6,3);
t119 = t31 + t106;
t118 = -t55 * mrSges(5,1) - t56 * mrSges(5,2);
t88 = t103 * pkin(3);
t117 = m(5) * t88;
t116 = -Ifges(6,5) + Ifges(5,4);
t115 = -t55 * t66 - t93 * t56;
t60 = pkin(3) * t66 + qJ(5);
t63 = -t93 * pkin(3) - pkin(4);
t114 = t60 * t55 - t63 * t56;
t89 = -m(5) / 0.2e1 - m(6) / 0.2e1;
t112 = t55 ^ 2;
t113 = -t56 ^ 2 - t112;
t57 = m(6) * t60 + mrSges(6,3);
t108 = m(6) / 0.2e1;
t107 = m(5) * pkin(3);
t105 = m(6) * t56;
t99 = t67 * mrSges(4,1);
t68 = -pkin(1) - pkin(6);
t58 = (-qJ(4) + t68) * t67;
t87 = t103 * t68;
t59 = -t103 * qJ(4) + t87;
t32 = t58 * t66 - t93 * t59;
t75 = t93 * t58 + t66 * t59;
t4 = (m(5) + m(6)) * (t32 * t55 + t56 * t75) + (mrSges(5,3) + mrSges(6,2)) * t113;
t96 = qJD(1) * t4;
t84 = t56 * mrSges(5,1) - t55 * mrSges(5,2);
t72 = -t103 * mrSges(4,2) - t31 - t84 - t99;
t9 = mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(5) * t64 + t106 - t72;
t95 = qJD(1) * t9;
t71 = (-t63 * t55 - t60 * t56) * t108 + (t93 * t55 - t56 * t66) * t107 / 0.2e1;
t28 = -t55 * pkin(4) + t56 * qJ(5) + t88;
t74 = t28 * t108 + t117 / 0.2e1;
t85 = -t55 * mrSges(6,1) + t56 * mrSges(6,3);
t5 = t118 - t71 + t74 + t85;
t94 = t5 * qJD(1);
t15 = t119 * t55;
t92 = qJD(1) * t15;
t80 = 0.2e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * t113;
t11 = t80 + t89;
t91 = t11 * qJD(1);
t27 = m(6) * t55;
t90 = t27 * qJD(1);
t1 = t84 * t88 + t25 * t85 + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t67) * t67 + ((Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3)) * t55 + t116 * t56) * t56 - t116 * t112 + (qJ(2) * mrSges(4,1) + (-Ifges(4,1) + Ifges(4,2)) * t67 - Ifges(4,4) * t103) * t103 + (t117 + t118) * t64 + t119 * t28;
t78 = t1 * qJD(1);
t76 = qJD(3) * t57;
t13 = -t56 * mrSges(6,2) + 0.2e1 * t75 * t108;
t10 = t80 - t89;
t8 = t71 + t74;
t2 = [qJD(2) * t9 + qJD(3) * t1 - qJD(4) * t4 + qJD(5) * t15, qJD(4) * t10 + t95, (-mrSges(4,2) * t87 - Ifges(4,5) * t67 - Ifges(4,6) * t103 - t68 * t99 + (-Ifges(6,4) - Ifges(5,5)) * t56 + (-Ifges(6,6) + Ifges(5,6)) * t55 - (t66 * t107 - mrSges(5,2) + t57) * t32 + (m(6) * t63 - t93 * t107 - mrSges(5,1) - mrSges(6,1)) * t75 - t115 * mrSges(5,3) * pkin(3) + t114 * mrSges(6,2)) * qJD(3) + t8 * qJD(4) + t13 * qJD(5) + t78, qJD(2) * t10 + qJD(3) * t8 - t96, qJD(3) * t13 + t92; qJD(4) * t11 - t95, 0, (-m(6) * t114 + t115 * t107 + t72) * qJD(3) + qJD(5) * t105, t91, qJD(3) * t105; -qJD(4) * t5 - t78, 0, t57 * qJD(5), -t94, t76; -qJD(2) * t11 + qJD(3) * t5 + qJD(5) * t27 + t96, -t91, t94, 0, t90; -qJD(4) * t27 - t92, 0, -t76, -t90, 0;];
Cq = t2;
