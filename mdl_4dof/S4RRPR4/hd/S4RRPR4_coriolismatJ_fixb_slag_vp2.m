% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:27
% EndTime: 2019-12-31 17:02:28
% DurationCPUTime: 0.43s
% Computational Cost: add. (1463->99), mult. (3061->133), div. (0->0), fcn. (2757->6), ass. (0->65)
t62 = sin(pkin(7));
t63 = cos(pkin(7));
t64 = sin(qJ(4));
t66 = cos(qJ(4));
t47 = -t64 * t62 + t66 * t63;
t98 = t47 / 0.2e1;
t48 = t66 * t62 + t64 * t63;
t97 = (Ifges(5,1) * t98 - Ifges(5,4) * t48 - Ifges(5,2) * t47 / 0.2e1) * t48;
t84 = t62 ^ 2 + t63 ^ 2;
t94 = mrSges(4,3) * t84;
t96 = -mrSges(3,2) + t94;
t76 = t84 * qJ(3);
t27 = t48 * mrSges(5,1) + t47 * mrSges(5,2);
t93 = (qJD(1) + qJD(2)) * t27;
t65 = sin(qJ(2));
t91 = pkin(2) * t65;
t90 = t65 * pkin(1);
t67 = cos(qJ(2));
t89 = t67 * pkin(1);
t56 = -t63 * pkin(3) - pkin(2);
t49 = t56 - t89;
t87 = t49 * t27;
t86 = t56 * t27;
t85 = Ifges(5,5) * t47 - Ifges(5,6) * t48;
t45 = Ifges(5,4) * t47;
t29 = -Ifges(5,2) * t48 + t45;
t32 = Ifges(5,1) * t48 + t45;
t75 = (t29 + t32) * t98 + t97;
t3 = t75 + t87;
t82 = t3 * qJD(1);
t55 = qJ(3) + t90;
t39 = (-pkin(6) - t55) * t62;
t59 = t63 * pkin(6);
t40 = t63 * t55 + t59;
t22 = t66 * t39 - t64 * t40;
t23 = t64 * t39 + t66 * t40;
t28 = -t47 * mrSges(5,1) + t48 * mrSges(5,2);
t37 = t48 * t89;
t38 = t47 * t89;
t51 = -t63 * mrSges(4,1) + t62 * mrSges(4,2);
t69 = (t37 * t48 + t38 * t47) * mrSges(5,3);
t77 = t84 * t55;
t5 = m(5) * (-t22 * t37 + t23 * t38 + t49 * t90) + t28 * t90 + m(4) * (-t91 + (t77 - t90) * t67) * pkin(1) + t51 * t90 - mrSges(3,1) * t90 + t69 + t96 * t89;
t81 = t5 * qJD(1);
t72 = (t47 ^ 2 + t48 ^ 2) * mrSges(5,3) + t94;
t10 = m(5) * (-t22 * t48 + t23 * t47) + m(4) * t77 + t72;
t80 = t10 * qJD(1);
t74 = (m(5) / 0.2e1 + m(4) / 0.2e1) * t90;
t70 = -t38 * mrSges(5,2) / 0.2e1 - t37 * mrSges(5,1) / 0.2e1;
t1 = -t97 + (-t32 / 0.2e1 - t29 / 0.2e1) * t47 + (-t49 / 0.2e1 - t56 / 0.2e1) * t27 + t70;
t4 = t75 + t86;
t73 = -t1 * qJD(1) + t4 * qJD(2);
t50 = (-pkin(6) - qJ(3)) * t62;
t52 = t63 * qJ(3) + t59;
t33 = t66 * t50 - t64 * t52;
t34 = t64 * t50 + t66 * t52;
t11 = m(5) * (-t33 * t48 + t34 * t47) + m(4) * t76 + t72;
t68 = -m(4) * (t77 + t76) / 0.2e1 - m(5) * ((-t22 - t33) * t48 + (t23 + t34) * t47) / 0.2e1 - t72;
t6 = t74 + t68;
t71 = -t6 * qJD(1) + t11 * qJD(2);
t26 = t27 * qJD(3);
t25 = t27 * qJD(4);
t7 = t74 - t68;
t2 = t87 / 0.2e1 + t86 / 0.2e1 + t70 + t75;
t8 = [t5 * qJD(2) + t10 * qJD(3) + t3 * qJD(4), t7 * qJD(3) + t2 * qJD(4) + t81 + (m(5) * (-t33 * t37 + t34 * t38) + t69 + (t96 * t67 + m(4) * (t67 * t76 - t91) + (m(5) * t56 - mrSges(3,1) + t28 + t51) * t65) * pkin(1)) * qJD(2), t7 * qJD(2) + t80, t82 + t2 * qJD(2) + (-t23 * mrSges(5,1) - t22 * mrSges(5,2) + t85) * qJD(4); -t6 * qJD(3) - t1 * qJD(4) - t81, t11 * qJD(3) + t4 * qJD(4), t71, (-t34 * mrSges(5,1) - t33 * mrSges(5,2) + t85) * qJD(4) + t73; t6 * qJD(2) + t25 - t80, t25 - t71, 0, t93; t1 * qJD(2) - t26 - t82, -t26 - t73, -t93, 0;];
Cq = t8;
