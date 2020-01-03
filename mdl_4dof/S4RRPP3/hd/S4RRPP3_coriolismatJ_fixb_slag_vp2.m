% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:24
% EndTime: 2019-12-31 16:57:25
% DurationCPUTime: 0.49s
% Computational Cost: add. (869->83), mult. (1835->104), div. (0->0), fcn. (1704->4), ass. (0->45)
t62 = sin(pkin(6));
t63 = cos(pkin(6));
t73 = sin(qJ(2));
t74 = cos(qJ(2));
t32 = t62 * t73 - t63 * t74;
t34 = -t62 * t74 - t63 * t73;
t45 = -t74 * pkin(2) - pkin(1);
t13 = t32 * pkin(3) + t34 * qJ(4) + t45;
t83 = m(5) * t13 + t32 * mrSges(5,1) + t34 * mrSges(5,3);
t82 = -t34 * mrSges(4,1) - t32 * mrSges(4,2);
t60 = t73 * pkin(2);
t81 = m(4) * t60;
t56 = t62 * pkin(2);
t39 = t56 + qJ(4);
t38 = m(5) * t39 + mrSges(5,3);
t78 = -Ifges(5,5) + Ifges(4,4);
t59 = t73 * pkin(5);
t52 = -t73 * qJ(3) - t59;
t48 = t74 * pkin(5);
t68 = t74 * qJ(3) + t48;
t21 = -t63 * t52 + t62 * t68;
t77 = t62 * t52 + t63 * t68;
t76 = m(5) / 0.2e1;
t75 = m(4) * pkin(2);
t72 = t32 * mrSges(5,2);
t16 = -t34 * pkin(3) + t32 * qJ(4) + t60;
t28 = t32 * mrSges(5,3);
t1 = -pkin(1) * (t73 * mrSges(3,1) + t74 * mrSges(3,2)) + t13 * t28 + (-t13 * mrSges(5,1) - mrSges(4,2) * t60 - t78 * t34) * t34 + (mrSges(4,1) * t60 + (Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,3)) * t34 + t78 * t32) * t32 + (-Ifges(3,2) + Ifges(3,1)) * t73 * t74 + (-t73 ^ 2 + t74 ^ 2) * Ifges(3,4) + (t81 + t82) * t45 + t83 * t16;
t67 = t1 * qJD(1);
t2 = (t32 ^ 2 + t34 ^ 2) * (mrSges(4,3) + mrSges(5,2)) + (m(4) + m(5)) * (-t21 * t34 - t32 * t77);
t66 = t2 * qJD(1);
t57 = t63 * pkin(2);
t44 = -t57 - pkin(3);
t49 = (-t39 * t32 - t44 * t34) * t76 + (-t62 * t32 + t63 * t34) * t75 / 0.2e1;
t51 = t16 * t76 + t81 / 0.2e1;
t3 = -t34 * mrSges(5,1) + t28 - t49 + t51 + t82;
t65 = t3 * qJD(1);
t5 = t83 * t34;
t64 = t5 * qJD(1);
t15 = m(5) * t34;
t61 = t15 * qJD(1);
t53 = t38 * qJD(2);
t6 = 0.2e1 * t77 * t76 - t72;
t4 = t49 + t51;
t7 = [t1 * qJD(2) + t2 * qJD(3) + t5 * qJD(4), t4 * qJD(3) + t6 * qJD(4) + t67 + (-mrSges(3,1) * t48 + mrSges(3,2) * t59 + Ifges(3,5) * t74 - Ifges(3,6) * t73 - t44 * t72 - (t62 * t75 - mrSges(4,2) + t38) * t21 + (m(5) * t44 - t63 * t75 - mrSges(4,1) - mrSges(5,1)) * t77 + (mrSges(4,3) * t57 - Ifges(5,4) - Ifges(4,5)) * t32 + (t39 * mrSges(5,2) + mrSges(4,3) * t56 + Ifges(4,6) - Ifges(5,6)) * t34) * qJD(2), t4 * qJD(2) + t66, t6 * qJD(2) + t64; -t3 * qJD(3) - t67, t38 * qJD(4), -t65, t53; t3 * qJD(2) + t15 * qJD(4) - t66, t65, 0, t61; -t15 * qJD(3) - t64, -t53, -t61, 0;];
Cq = t7;
