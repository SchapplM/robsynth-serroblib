% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:22
% EndTime: 2019-12-31 17:53:23
% DurationCPUTime: 0.36s
% Computational Cost: add. (1040->82), mult. (2215->107), div. (0->0), fcn. (2052->4), ass. (0->45)
t63 = cos(pkin(7));
t62 = sin(pkin(7));
t73 = -t63 * pkin(2) - t62 * qJ(3) - pkin(1);
t45 = t63 * pkin(3) - t73;
t64 = cos(qJ(4));
t86 = sin(qJ(4));
t49 = t62 * t86 + t63 * t64;
t50 = t62 * t64 - t63 * t86;
t15 = t49 * pkin(4) - t50 * qJ(5) + t45;
t98 = m(6) * t15 + t49 * mrSges(6,1) - t50 * mrSges(6,3);
t97 = m(6) + m(5);
t95 = mrSges(5,1) + mrSges(6,1);
t94 = -Ifges(6,5) + Ifges(5,4);
t93 = m(6) / 0.2e1 + m(5) / 0.2e1;
t57 = m(6) * qJ(5) + mrSges(6,3);
t92 = -m(6) * pkin(4) - t95;
t91 = -mrSges(5,2) + t57;
t87 = t50 * pkin(4);
t84 = -pkin(6) + qJ(2);
t82 = qJ(5) * t49;
t66 = t82 + t87;
t1 = (t45 * mrSges(5,1) + t15 * mrSges(6,1) - t94 * t50) * t50 + (-t45 * mrSges(5,2) + t15 * mrSges(6,3) + (-Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,3)) * t50 + t94 * t49) * t49 + t98 * t66;
t81 = t1 * qJD(1);
t52 = t84 * t63;
t70 = t84 * t62;
t26 = t52 * t86 - t64 * t70;
t27 = t64 * t52 + t70 * t86;
t2 = (t49 ^ 2 + t50 ^ 2) * (-mrSges(5,3) - mrSges(6,2)) + t97 * (-t26 * t50 + t27 * t49) + (mrSges(3,3) + mrSges(4,2) + 0.4e1 * (m(3) / 0.4e1 + m(4) / 0.4e1) * qJ(2)) * (t62 ^ 2 + t63 ^ 2);
t80 = t2 * qJD(1);
t5 = t95 * t50 + (-mrSges(5,2) + mrSges(6,3)) * t49 + 0.2e1 * (t87 / 0.4e1 + t82 / 0.4e1 + t66 / 0.4e1) * m(6);
t79 = t5 * qJD(1);
t7 = (-m(4) * t73 + m(5) * t45 + t63 * mrSges(4,1) + t49 * mrSges(5,1) + t50 * mrSges(5,2) + t62 * mrSges(4,3) + t98) * t62;
t78 = t7 * qJD(1);
t8 = t98 * t50;
t77 = t8 * qJD(1);
t65 = t93 * (t86 * t49 + t64 * t50);
t9 = (m(4) + t93) * t62 + t65;
t76 = t9 * qJD(1);
t20 = m(6) * t50;
t75 = t20 * qJD(1);
t74 = t57 * qJD(4);
t72 = m(6) * t86;
t13 = m(6) * t27 - t49 * mrSges(6,2);
t10 = t65 - t97 * t62 / 0.2e1;
t3 = [t2 * qJD(2) + t7 * qJD(3) + t1 * qJD(4) - t8 * qJD(5), t10 * qJD(3) + t80, t10 * qJD(2) + t78, t13 * qJD(5) + t81 + (-t91 * t26 + t92 * t27 + (-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t50 + (pkin(4) * mrSges(6,2) - Ifges(6,4) - Ifges(5,5)) * t49) * qJD(4), t13 * qJD(4) - t77; -t9 * qJD(3) - t5 * qJD(4) + t20 * qJD(5) - t80, 0, -t76, -t79, t75; t9 * qJD(2) - t78, t76, 0, (t91 * t64 + t92 * t86) * qJD(4) + qJD(5) * t72, qJD(4) * t72; t5 * qJD(2) - t81, t79, 0, t57 * qJD(5), t74; -t20 * qJD(2) + t77, -t75, 0, -t74, 0;];
Cq = t3;
