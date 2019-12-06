% Calculate vector of cutting forces with Newton-Euler
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:42
% EndTime: 2019-12-05 15:32:43
% DurationCPUTime: 0.46s
% Computational Cost: add. (3034->102), mult. (5280->130), div. (0->0), fcn. (2789->8), ass. (0->54)
t57 = sin(qJ(4));
t59 = cos(qJ(4));
t73 = qJD(2) * qJD(4);
t66 = t59 * t73;
t35 = t57 * qJDD(2) + t66;
t36 = t59 * qJDD(2) - t57 * t73;
t75 = qJD(2) * t57;
t43 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t75;
t74 = qJD(2) * t59;
t44 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t74;
t45 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t74;
t61 = qJD(2) ^ 2;
t54 = sin(pkin(7));
t56 = cos(pkin(7));
t40 = -t56 * g(1) - t54 * g(2);
t52 = -g(3) + qJDD(1);
t58 = sin(qJ(2));
t60 = cos(qJ(2));
t64 = -t58 * t40 + t60 * t52;
t23 = qJDD(2) * pkin(2) + t64;
t77 = t60 * t40 + t58 * t52;
t24 = -t61 * pkin(2) + t77;
t53 = sin(pkin(8));
t55 = cos(pkin(8));
t65 = t55 * t23 - t53 * t24;
t63 = -qJDD(2) * pkin(3) - t65;
t41 = qJD(4) * pkin(4) - qJ(5) * t75;
t42 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t75;
t51 = t59 ^ 2;
t67 = m(6) * (t41 * t75 - t36 * pkin(4) + qJDD(5) + (-qJ(5) * t51 - pkin(6)) * t61 + t63) + t42 * t75 + t35 * mrSges(6,2);
t85 = -(-t57 * t43 + (t44 + t45) * t59) * qJD(2) - (mrSges(5,1) + mrSges(6,1)) * t36 + m(5) * (-t61 * pkin(6) + t63) + t35 * mrSges(5,2) + t67;
t82 = pkin(4) * t61;
t78 = t53 * t23 + t55 * t24;
t19 = -t61 * pkin(3) + qJDD(2) * pkin(6) + t78;
t39 = t54 * g(1) - t56 * g(2);
t38 = qJDD(3) - t39;
t79 = t59 * t19 + t57 * t38;
t72 = qJD(2) * qJD(5);
t27 = t59 * t38;
t33 = (-mrSges(6,1) * t59 + mrSges(6,2) * t57) * qJD(2);
t34 = (-mrSges(5,1) * t59 + mrSges(5,2) * t57) * qJD(2);
t69 = m(6) * (qJDD(4) * pkin(4) + t27 + (-t35 + t66) * qJ(5) + (t59 * t82 - t19 - 0.2e1 * t72) * t57) + qJD(4) * t44 + qJDD(4) * mrSges(6,1);
t11 = m(5) * (-t57 * t19 + t27) + qJDD(4) * mrSges(5,1) + qJD(4) * t45 + (-mrSges(5,3) - mrSges(6,3)) * t35 + (-t33 - t34) * t75 + t69;
t68 = m(6) * (t36 * qJ(5) - qJD(4) * t41 - t51 * t82 + 0.2e1 * t59 * t72 + t79) + t33 * t74 + t36 * mrSges(6,3);
t12 = m(5) * t79 + t36 * mrSges(5,3) + t34 * t74 + (-mrSges(5,2) - mrSges(6,2)) * qJDD(4) + (-t43 - t42) * qJD(4) + t68;
t6 = m(4) * t78 - t61 * mrSges(4,1) - qJDD(2) * mrSges(4,2) - t57 * t11 + t59 * t12;
t8 = m(4) * t65 + qJDD(2) * mrSges(4,1) - t61 * mrSges(4,2) - t85;
t4 = m(3) * t64 + qJDD(2) * mrSges(3,1) - t61 * mrSges(3,2) + t53 * t6 + t55 * t8;
t5 = m(3) * t77 - t61 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t53 * t8 + t55 * t6;
t71 = m(2) * t52 + t60 * t4 + t58 * t5;
t70 = m(4) * t38 + t59 * t11 + t57 * t12;
t7 = (m(2) + m(3)) * t39 - t70;
t1 = m(2) * t40 - t58 * t4 + t60 * t5;
t2 = [-m(1) * g(1) + t56 * t1 - t54 * t7, t1, t5, t6, t12, -qJDD(4) * mrSges(6,2) - qJD(4) * t42 + t68; -m(1) * g(2) + t54 * t1 + t56 * t7, t7, t4, t8, t11, -t35 * mrSges(6,3) - t33 * t75 + t69; -m(1) * g(3) + t71, t71, -m(3) * t39 + t70, t70, t85, -t36 * mrSges(6,1) - t44 * t74 + t67;];
f_new = t2;
