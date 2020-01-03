% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:42
% EndTime: 2019-12-31 19:27:43
% DurationCPUTime: 0.52s
% Computational Cost: add. (1076->117), mult. (2028->157), div. (0->0), fcn. (1385->6), ass. (0->78)
t124 = -Ifges(6,1) + Ifges(6,2);
t72 = cos(qJ(5));
t115 = -t72 / 0.2e1;
t67 = sin(pkin(8));
t70 = sin(qJ(5));
t87 = t72 * mrSges(6,1) - t70 * mrSges(6,2);
t123 = t67 * t87;
t68 = cos(pkin(8));
t103 = t72 * mrSges(6,2);
t106 = t70 * mrSges(6,1);
t86 = -t103 - t106;
t27 = t68 * t86;
t122 = (qJD(1) + qJD(2)) * t27;
t66 = t72 ^ 2;
t99 = t70 ^ 2 + t66;
t71 = sin(qJ(2));
t73 = cos(qJ(2));
t120 = -t73 * mrSges(3,2) + (-mrSges(3,1) - mrSges(4,1)) * t71;
t119 = mrSges(6,3) * t99;
t118 = m(5) / 0.2e1;
t117 = m(6) / 0.2e1;
t38 = (t67 * t71 + t68 * t73) * pkin(1);
t116 = -t38 / 0.2e1;
t114 = t71 * pkin(1);
t113 = t73 * pkin(1);
t112 = Ifges(6,4) * t70;
t111 = Ifges(6,4) * t72;
t93 = -pkin(2) - t113;
t56 = -pkin(3) + t93;
t50 = t67 * t56;
t58 = qJ(3) + t114;
t109 = t67 * t58;
t37 = t67 * t113 - t68 * t114;
t107 = t68 * t37;
t81 = -t106 / 0.2e1 - t103 / 0.2e1;
t28 = -t27 / 0.2e1 + t81 * t68;
t101 = t28 * qJD(5) + m(6) * (-0.1e1 + t99) * t68 * t67 * qJD(3);
t100 = t68 * t58 + t50;
t74 = -pkin(2) - pkin(3);
t48 = t68 * qJ(3) + t67 * t74;
t33 = -pkin(7) + t100;
t77 = -mrSges(4,3) * t113 + (-mrSges(5,2) + t119) * t38 + (-mrSges(5,1) - t87) * t37;
t85 = t68 * t56 - t109;
t82 = pkin(4) - t85;
t92 = t99 * t38;
t3 = -m(6) * (t33 * t92 + t82 * t37) - m(5) * (t100 * t38 - t85 * t37) + t77 + (-m(4) * (t58 * t73 + t93 * t71) - t120) * pkin(1);
t98 = t3 * qJD(1);
t88 = -t67 * mrSges(5,1) - t68 * mrSges(5,2) - mrSges(4,3) - t123;
t79 = t68 * t119 + t88;
t90 = t99 * t68;
t4 = -m(6) * (t33 * t90 + t82 * t67) - m(5) * (t100 * t68 - t85 * t67) - m(4) * t58 + t79;
t97 = t4 * qJD(1);
t26 = t82 * t86;
t80 = t124 * t72 + t112;
t78 = -t66 * Ifges(6,4) + t80 * t70;
t11 = -t26 + t78;
t96 = t11 * qJD(1);
t47 = -t67 * qJ(3) + t68 * t74;
t45 = pkin(4) - t47;
t32 = t45 * t86;
t94 = -t26 / 0.2e1 - t32 / 0.2e1;
t46 = -pkin(7) + t48;
t91 = t99 * t46;
t89 = -Ifges(6,5) * t72 + Ifges(6,6) * t70;
t75 = m(4) * (0.4e1 * qJ(3) + 0.2e1 * t114) / 0.4e1 + ((t48 + t100) * t68 + (-t47 - t85) * t67) * t118 + ((pkin(4) + t45 + t109) * t67 + (t99 * t33 - t50 + t91) * t68) * t117 - t79;
t76 = -m(5) * (t67 * t38 - t107) / 0.2e1 - m(6) * (t67 * t92 - t107) / 0.2e1 - m(4) * t114 / 0.2e1;
t2 = t75 + t76;
t8 = -m(4) * qJ(3) + mrSges(6,3) * t90 - m(6) * (t45 * t67 + t46 * t90) - m(5) * (-t47 * t67 + t48 * t68) + t88;
t84 = -t2 * qJD(1) + t8 * qJD(2);
t14 = -t32 + t78;
t5 = (mrSges(6,2) * t116 - t111) * t72 + (mrSges(6,1) * t116 + t80) * t70 + t94;
t83 = -t5 * qJD(1) - t14 * qJD(2);
t24 = t28 * qJD(3);
t22 = t27 * qJD(5);
t21 = t27 * qJD(3);
t6 = -0.2e1 * t111 * t115 + t81 * t38 - t94 + (0.2e1 * t124 * t115 - t112) * t70;
t1 = t75 - t76;
t7 = [-t3 * qJD(2) - t4 * qJD(3) - t11 * qJD(5), t1 * qJD(3) + t6 * qJD(5) - t98 + (-t77 + 0.2e1 * (t45 * t37 + t38 * t91) * t117 + 0.2e1 * (-t47 * t37 + t48 * t38) * t118 + (m(4) * (-pkin(2) * t71 + qJ(3) * t73) + t120) * pkin(1)) * qJD(2), t1 * qJD(2) + t101 - t97, 0, -t96 + t6 * qJD(2) + t24 + (-t87 * t33 + t89) * qJD(5); t2 * qJD(3) - t5 * qJD(5) + t98, -t8 * qJD(3) - t14 * qJD(5), -t84 + t101, 0, t24 + (-t87 * t46 + t89) * qJD(5) + t83; -t2 * qJD(2) - t22 + t97, -t22 + t84, 0, 0, -qJD(5) * t123 - t122; 0, 0, 0, 0, t86 * qJD(5); t5 * qJD(2) + t21 + t96, t21 - t83, t122, 0, 0;];
Cq = t7;
