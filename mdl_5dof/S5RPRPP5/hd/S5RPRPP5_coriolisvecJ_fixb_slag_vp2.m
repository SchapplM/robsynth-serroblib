% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:08
% EndTime: 2019-12-31 18:16:10
% DurationCPUTime: 0.96s
% Computational Cost: add. (593->193), mult. (1312->240), div. (0->0), fcn. (448->2), ass. (0->77)
t97 = Ifges(6,4) + Ifges(5,5);
t95 = ((qJ(2) * (m(3) + m(4)) + mrSges(3,3)) * qJD(1));
t50 = -pkin(1) - pkin(6);
t41 = qJD(1) * t50 + qJD(2);
t48 = cos(qJ(3));
t68 = qJ(5) * qJD(1);
t13 = (t41 + t68) * t48;
t94 = -t13 + qJD(4);
t93 = (-mrSges(4,1) - mrSges(5,1)) * qJD(3);
t90 = -qJD(3) / 0.2e1;
t89 = qJD(3) / 0.2e1;
t69 = t48 * qJD(1);
t88 = -t97 * t69 / 0.2e1;
t49 = -pkin(3) - pkin(4);
t85 = Ifges(4,4) * t48;
t84 = t41 * t48;
t47 = sin(qJ(3));
t31 = t47 * t41;
t83 = -mrSges(5,2) + mrSges(6,3);
t73 = qJD(3) * t41;
t15 = qJD(3) * qJD(4) + t48 * t73;
t28 = (mrSges(5,1) * t47 - mrSges(5,3) * t48) * qJD(1);
t29 = (-mrSges(6,1) * t47 + t48 * mrSges(6,2)) * qJD(1);
t82 = t28 - t29;
t24 = qJD(3) * qJ(4) + t31;
t70 = t47 * qJD(1);
t35 = qJD(3) * mrSges(6,2) + mrSges(6,3) * t70;
t40 = -mrSges(5,2) * t70 + qJD(3) * mrSges(5,3);
t81 = t35 + t40;
t74 = qJD(3) * mrSges(4,2);
t80 = -mrSges(4,3) * t70 + t40 - t74;
t79 = (mrSges(5,2) + mrSges(4,3)) * t69 + t93;
t78 = (qJ(2) * mrSges(4,1));
t77 = qJ(2) * mrSges(4,2);
t76 = qJ(4) * t47;
t75 = qJ(5) + t50;
t72 = qJD(5) * t47;
t71 = qJD(5) * t48;
t67 = qJ(5) * qJD(3);
t66 = qJD(1) * qJD(2);
t65 = t47 * t73;
t16 = -qJD(3) * pkin(3) + qJD(4) - t84;
t64 = -t16 + t84;
t63 = qJ(4) * t48 - qJ(2);
t62 = qJD(3) * t75;
t60 = -qJD(4) * t48 + qJD(2);
t59 = 0.3e1 / 0.2e1 * Ifges(6,4) + 0.3e1 / 0.2e1 * Ifges(5,5) - 0.3e1 / 0.2e1 * Ifges(4,4);
t58 = -Ifges(4,5) / 0.2e1 + Ifges(6,5) / 0.2e1 - Ifges(5,4) / 0.2e1;
t57 = -Ifges(4,6) / 0.2e1 - Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t56 = pkin(3) * t48 + t76;
t55 = t16 * t47 + t24 * t48;
t34 = pkin(3) * t47 - t63;
t54 = t48 * t49 - t76;
t25 = t47 * t49 + t63;
t11 = qJD(3) * t56 + t60;
t17 = t34 * qJD(1);
t4 = qJD(3) * t49 + t94;
t7 = qJD(1) * t25 + qJD(5);
t53 = t7 * mrSges(6,2) - t17 * mrSges(5,3) - t4 * mrSges(6,3) + Ifges(6,5) * t90 + (Ifges(5,4) + Ifges(4,5)) * t89 + ((Ifges(4,1) + Ifges(5,1) + Ifges(6,1)) * t48 + (-Ifges(4,4) + t97) * t47) * qJD(1) / 0.2e1;
t6 = qJD(3) * t54 - t60;
t42 = t47 * t68;
t8 = t42 + t24;
t52 = t17 * mrSges(5,1) + t8 * mrSges(6,3) + Ifges(5,6) * t89 - (-Ifges(4,2) * t47 + t85) * qJD(1) / 0.2e1 - t24 * mrSges(5,2) - t7 * mrSges(6,1) - t88 + (Ifges(5,3) + Ifges(6,2)) * t70 / 0.2e1 + (Ifges(6,6) + Ifges(4,6)) * t90;
t37 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t69;
t33 = t75 * t48;
t32 = t75 * t47;
t30 = qJD(1) * (mrSges(4,1) * t47 + mrSges(4,2) * t48);
t27 = t56 * qJD(1);
t14 = t54 * qJD(1);
t12 = t31 + t42;
t10 = t48 * t62 + t72;
t9 = t47 * t62 - t71;
t5 = t11 * qJD(1);
t3 = t65 + (t47 * t67 - t71) * qJD(1);
t2 = (t48 * t67 + t72) * qJD(1) + t15;
t1 = t6 * qJD(1);
t18 = [t10 * t35 + t11 * t28 + t6 * t29 + t9 * t37 + m(5) * (t11 * t17 + t34 * t5) + m(6) * (t1 * t25 + t10 * t8 + t2 * t32 - t3 * t33 + t4 * t9 + t6 * t7) + (t30 + (2 * t95)) * qJD(2) + (mrSges(4,2) * t66 + t1 * mrSges(6,2) - t5 * mrSges(5,3) - t3 * mrSges(6,3) + (t57 * qJD(3) + (m(5) * t24 + t80) * t50 + (t34 * mrSges(5,1) - t25 * mrSges(6,1) + t32 * mrSges(6,3) + t48 * t59 + (2 * t78)) * qJD(1) + t52) * qJD(3)) * t48 + (mrSges(4,1) * t66 + t5 * mrSges(5,1) - t1 * mrSges(6,1) + t2 * mrSges(6,3) + (m(5) * t50 - mrSges(5,2)) * t15 + (t58 * qJD(3) + t64 * mrSges(5,2) + (-m(5) * t64 + t79) * t50 + (-t33 * mrSges(6,3) - 0.2e1 * t77 + t34 * mrSges(5,3) - t25 * mrSges(6,2) - t59 * t47 + (-0.3e1 / 0.2e1 * Ifges(6,1) - 0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(6,2)) * t48) * qJD(1) - t53) * qJD(3)) * t47; m(5) * t15 * t47 + m(6) * (t2 * t47 - t3 * t48) + ((t35 + t80) * t48 + (t37 + t79) * t47 + m(5) * (-t31 * t48 + t55) + m(6) * (t4 * t47 + t48 * t8)) * qJD(3) + (-m(5) * t17 + m(6) * t7 - t30 - t82 - t95) * qJD(1); -t3 * mrSges(6,1) + t2 * mrSges(6,2) + t15 * mrSges(5,3) - t12 * t37 - t13 * t35 - t14 * t29 - t27 * t28 + t81 * qJD(4) + ((-t74 - t80) * t48 + (t93 - t79) * t47) * t41 + (((-t78 + t85 / 0.2e1) * qJD(1) + (qJ(4) * t83 + t57) * qJD(3) - t52 + t88) * t48 + (t16 * mrSges(5,2) + (t77 + (-Ifges(4,4) / 0.2e1 + Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t47) * qJD(1) + (Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1 - Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t69 + (mrSges(5,2) * pkin(3) + mrSges(6,3) * t49 + t58) * qJD(3) + t53) * t47) * qJD(1) + (t2 * qJ(4) - t4 * t12 - t7 * t14 + t3 * t49 + t8 * t94) * m(6) + (-pkin(3) * t65 + qJ(4) * t15 + qJD(4) * t24 - t17 * t27 - t41 * t55) * m(5); t82 * t69 + ((m(5) * t41 + qJD(1) * t83) * t47 - t81) * qJD(3) - m(5) * (qJD(3) * t24 - t17 * t69) + (-t8 * qJD(3) - t69 * t7 + t3) * m(6); m(6) * t1 + (t48 * t37 - t47 * t35 - m(6) * (-t4 * t48 + t47 * t8) + (-mrSges(6,1) * t48 - mrSges(6,2) * t47) * qJD(3)) * qJD(1);];
tauc = t18(:);
