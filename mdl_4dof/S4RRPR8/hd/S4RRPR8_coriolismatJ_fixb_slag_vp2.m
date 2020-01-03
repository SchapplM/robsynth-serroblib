% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:49
% EndTime: 2019-12-31 17:07:51
% DurationCPUTime: 0.54s
% Computational Cost: add. (1206->83), mult. (2322->106), div. (0->0), fcn. (1980->4), ass. (0->46)
t61 = sin(qJ(2));
t62 = cos(qJ(4));
t63 = cos(qJ(2));
t91 = sin(qJ(4));
t46 = -t61 * t91 - t62 * t63;
t94 = t46 / 0.2e1;
t112 = 0.2e1 * t94;
t102 = pkin(5) - pkin(6);
t53 = t102 * t61;
t54 = t102 * t63;
t29 = t53 * t62 - t54 * t91;
t47 = t61 * t62 - t63 * t91;
t67 = t53 * t91 + t54 * t62;
t3 = t67 * mrSges(5,1) + t29 * mrSges(5,2) - Ifges(5,5) * t46 + Ifges(5,6) * t47;
t109 = t3 * qJD(4);
t82 = t61 * qJ(3);
t92 = pkin(2) + pkin(3);
t40 = t63 * t92 + pkin(1) + t82;
t108 = m(5) * t40 - mrSges(5,1) * t46 + mrSges(5,2) * t47;
t66 = m(5) * (-t29 * t91 + t62 * t67);
t65 = t66 / 0.2e1;
t68 = mrSges(5,1) * t91 + mrSges(5,2) * t62;
t104 = t68 * qJD(4);
t103 = (-t40 * mrSges(5,2) - 0.2e1 * Ifges(5,4) * t94) * t46 + (-t40 * mrSges(5,1) + Ifges(5,4) * t47 + (-Ifges(5,1) + Ifges(5,2)) * t112) * t47;
t99 = Ifges(3,4) - Ifges(4,5);
t72 = -pkin(2) * t63 - t82;
t50 = -pkin(1) + t72;
t51 = -t63 * mrSges(4,1) - t61 * mrSges(4,3);
t75 = m(4) * t50 + t51;
t81 = t63 * qJ(3);
t1 = t75 * (pkin(2) * t61 - t81) + (-pkin(1) * mrSges(3,2) - t50 * mrSges(4,3) + t99 * t63) * t63 + (-pkin(1) * mrSges(3,1) + t50 * mrSges(4,1) + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t63 - t99 * t61) * t61 + t103 + t108 * (-t61 * t92 + t81);
t84 = t1 * qJD(1);
t83 = t103 * qJD(1);
t10 = (-t75 + t108) * t61;
t80 = t10 * qJD(1);
t48 = -qJ(3) * t91 - t62 * t92;
t49 = qJ(3) * t62 - t91 * t92;
t11 = mrSges(5,1) * t49 + mrSges(5,2) * t48;
t79 = t11 * qJD(4);
t78 = t68 * qJD(2);
t71 = t11 * qJD(2);
t15 = mrSges(4,3) + m(4) * qJ(3) + m(5) * (-t48 * t91 + t49 * t62) + t68;
t6 = t65 - t66 / 0.2e1;
t70 = qJD(1) * t6 + qJD(2) * t15;
t5 = 0.2e1 * t65 + (m(4) * pkin(5) + mrSges(4,2)) * t63 + (t112 * t62 + t47 * t91) * mrSges(5,3);
t2 = [qJD(2) * t1 + qJD(3) * t10 - qJD(4) * t103, t5 * qJD(3) + t109 + t84 + (m(5) * (-t29 * t49 + t48 * t67) + (t46 * t48 + t47 * t49) * mrSges(5,3) + (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t63 + (-qJ(3) * mrSges(4,2) - Ifges(3,6) + Ifges(4,6)) * t61 + (m(4) * t72 - t63 * mrSges(3,1) + t61 * mrSges(3,2) + t51) * pkin(5) - t3) * qJD(2), qJD(2) * t5 + t80, t3 * qJD(2) - t109 - t83; qJD(3) * t6 - t84, qJD(3) * t15 + t79, t70, t71 - t79; -qJD(2) * t6 - t80, -t70 + t104, 0, t78 - t104; t83, -qJD(3) * t68 - t71, -t78, 0;];
Cq = t2;
