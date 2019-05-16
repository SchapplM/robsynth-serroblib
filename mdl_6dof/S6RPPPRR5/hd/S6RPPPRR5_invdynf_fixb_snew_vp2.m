% Calculate vector of cutting forces with Newton-Euler
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPPRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:50:08
% EndTime: 2019-05-05 13:50:10
% DurationCPUTime: 0.66s
% Computational Cost: add. (6549->130), mult. (11040->153), div. (0->0), fcn. (4520->8), ass. (0->71)
t97 = -2 * qJD(1);
t58 = sin(qJ(1));
t61 = cos(qJ(1));
t73 = -t61 * g(1) - t58 * g(2);
t96 = -qJDD(1) * qJ(2) + (qJD(2) * t97) - t73;
t63 = qJD(1) ^ 2;
t79 = t58 * g(1) - t61 * g(2);
t72 = qJDD(2) - t79;
t90 = -pkin(1) - qJ(3);
t67 = (qJD(3) * t97) + t90 * qJDD(1) + t72;
t26 = (-pkin(3) - qJ(2)) * t63 + t67;
t66 = t90 * t63 + qJDD(3) - t96;
t27 = qJDD(1) * pkin(3) + t66;
t54 = sin(pkin(9));
t55 = cos(pkin(9));
t75 = -t54 * t26 + t55 * t27;
t18 = -qJDD(1) * pkin(4) - (t63 * pkin(7)) - t75;
t57 = sin(qJ(5));
t60 = cos(qJ(5));
t83 = qJD(1) * qJD(5);
t77 = t60 * t83;
t41 = t57 * qJDD(1) + t77;
t78 = t57 * t83;
t42 = t60 * qJDD(1) - t78;
t14 = (-t41 - t77) * pkin(8) + (-t42 + t78) * pkin(5) + t18;
t40 = (-pkin(5) * t60 - pkin(8) * t57) * qJD(1);
t62 = qJD(5) ^ 2;
t85 = t60 * qJD(1);
t88 = t55 * t26 + t54 * t27;
t19 = -t63 * pkin(4) + qJDD(1) * pkin(7) + t88;
t53 = -g(3) + qJDD(4);
t89 = t60 * t19 + t57 * t53;
t16 = -t62 * pkin(5) + qJDD(5) * pkin(8) + t40 * t85 + t89;
t56 = sin(qJ(6));
t59 = cos(qJ(6));
t86 = qJD(1) * t57;
t37 = t59 * qJD(5) - t56 * t86;
t21 = t37 * qJD(6) + t56 * qJDD(5) + t59 * t41;
t38 = t56 * qJD(5) + t59 * t86;
t25 = -t37 * mrSges(7,1) + t38 * mrSges(7,2);
t45 = qJD(6) - t85;
t32 = -t45 * mrSges(7,2) + t37 * mrSges(7,3);
t36 = qJDD(6) - t42;
t12 = m(7) * (t59 * t14 - t56 * t16) - t21 * mrSges(7,3) + t36 * mrSges(7,1) - t38 * t25 + t45 * t32;
t20 = -t38 * qJD(6) + t59 * qJDD(5) - t56 * t41;
t33 = t45 * mrSges(7,1) - t38 * mrSges(7,3);
t13 = m(7) * (t56 * t14 + t59 * t16) + t20 * mrSges(7,3) - t36 * mrSges(7,2) + t37 * t25 - t45 * t33;
t43 = (qJD(5) * mrSges(6,1)) - mrSges(6,3) * t86;
t44 = -(qJD(5) * mrSges(6,2)) + mrSges(6,3) * t85;
t95 = m(6) * t18 - t42 * mrSges(6,1) + t41 * mrSges(6,2) + t59 * t12 + t56 * t13 + (t57 * t43 - t60 * t44) * qJD(1);
t94 = -m(3) - m(4);
t93 = t60 * t53;
t92 = -mrSges(2,2) + mrSges(3,3);
t91 = (mrSges(3,2) - mrSges(4,3));
t87 = t63 * qJ(2);
t84 = -m(2) + t94;
t39 = (-mrSges(6,1) * t60 + mrSges(6,2) * t57) * qJD(1);
t64 = m(7) * (-qJDD(5) * pkin(5) - t62 * pkin(8) - t93 + (qJD(1) * t40 + t19) * t57) - t20 * mrSges(7,1) + t21 * mrSges(7,2) - t37 * t32 + t38 * t33;
t11 = m(6) * (-t57 * t19 + t93) - t41 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t39 * t86 + qJD(5) * t44 - t64;
t9 = m(6) * t89 - qJDD(5) * mrSges(6,2) + t42 * mrSges(6,3) - qJD(5) * t43 - t56 * t12 + t59 * t13 + t39 * t85;
t81 = m(5) * t53 + t60 * t11 + t57 * t9;
t80 = mrSges(2,1) - t91;
t5 = m(5) * t88 - (t63 * mrSges(5,1)) - qJDD(1) * mrSges(5,2) - t57 * t11 + t60 * t9;
t7 = m(5) * t75 + qJDD(1) * mrSges(5,1) - t63 * mrSges(5,2) - t95;
t76 = m(4) * t66 + qJDD(1) * mrSges(4,1) + t54 * t5 + t55 * t7;
t74 = -t55 * t5 + t54 * t7 - m(4) * (t67 - t87);
t69 = m(3) * ((t63 * pkin(1)) + t96) - t76;
t68 = m(3) * (-qJDD(1) * pkin(1) + t72 - t87) - t74;
t2 = m(2) * t79 + (mrSges(4,1) + t92) * t63 + t80 * qJDD(1) - t68;
t1 = m(2) * t73 + t92 * qJDD(1) - t80 * t63 - t69;
t3 = [-m(1) * g(1) + t61 * t1 - t58 * t2, t1, t94 * g(3) + t81, -t63 * mrSges(4,1) - qJDD(1) * mrSges(4,3) - t74, t5, t9, t13; -m(1) * g(2) + t58 * t1 + t61 * t2, t2, -qJDD(1) * mrSges(3,3) - (t91 * t63) + t69, m(4) * g(3) - t81, t7, t11, t12; (-m(1) + t84) * g(3) + t81, t84 * g(3) + t81, (-mrSges(4,1) - mrSges(3,3)) * t63 + t91 * qJDD(1) + t68, -t63 * mrSges(4,3) + t76, t81, t95, t64;];
f_new  = t3;
