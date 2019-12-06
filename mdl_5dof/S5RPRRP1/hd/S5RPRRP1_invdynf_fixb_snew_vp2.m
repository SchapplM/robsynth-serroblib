% Calculate vector of cutting forces with Newton-Euler
% S5RPRRP1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:16
% EndTime: 2019-12-05 17:59:19
% DurationCPUTime: 0.61s
% Computational Cost: add. (4389->135), mult. (8724->165), div. (0->0), fcn. (5015->6), ass. (0->63)
t69 = qJD(1) ^ 2;
t64 = sin(qJ(3));
t67 = cos(qJ(3));
t88 = qJD(1) * qJD(3);
t50 = -t64 * qJDD(1) - t67 * t88;
t80 = t64 * t88;
t51 = t67 * qJDD(1) - t80;
t90 = qJD(1) * t64;
t52 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t90;
t89 = qJD(1) * t67;
t53 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t89;
t63 = sin(qJ(4));
t66 = cos(qJ(4));
t44 = (-t63 * t64 + t66 * t67) * qJD(1);
t23 = -t44 * qJD(4) + t66 * t50 - t63 * t51;
t43 = (-t63 * t67 - t64 * t66) * qJD(1);
t24 = t43 * qJD(4) + t63 * t50 + t66 * t51;
t60 = qJD(3) + qJD(4);
t33 = -t60 * mrSges(6,2) + t43 * mrSges(6,3);
t34 = -t60 * mrSges(5,2) + t43 * mrSges(5,3);
t37 = t60 * mrSges(5,1) - t44 * mrSges(5,3);
t54 = qJD(3) * pkin(3) - pkin(7) * t89;
t62 = t64 ^ 2;
t65 = sin(qJ(1));
t68 = cos(qJ(1));
t78 = -t68 * g(1) - t65 * g(2);
t75 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t78;
t99 = -pkin(1) - pkin(6);
t71 = -t50 * pkin(3) + t54 * t89 + (-pkin(7) * t62 + t99) * t69 + t75;
t35 = t60 * pkin(4) - t44 * qJ(5);
t36 = t60 * mrSges(6,1) - t44 * mrSges(6,3);
t42 = t43 ^ 2;
t84 = m(6) * (-t23 * pkin(4) - t42 * qJ(5) + t44 * t35 + qJDD(5) + t71) + t24 * mrSges(6,2) + t44 * t36;
t72 = m(5) * t71 + t24 * mrSges(5,2) - (mrSges(5,1) + mrSges(6,1)) * t23 + t44 * t37 + (-t34 - t33) * t43 + t84;
t70 = m(4) * (t99 * t69 + t75) + t52 * t90 + t53 * t89 + t51 * mrSges(4,2) + t72 - t50 * mrSges(4,1);
t101 = m(3) * (t69 * pkin(1) - t75) - t70;
t100 = -m(2) - m(3);
t96 = mrSges(2,1) - mrSges(3,2);
t94 = -mrSges(2,2) + mrSges(3,3);
t82 = t65 * g(1) - t68 * g(2);
t73 = -t69 * qJ(2) + qJDD(2) - t82;
t39 = t99 * qJDD(1) + t73;
t91 = t64 * g(3) + t67 * t39;
t15 = (-t51 - t80) * pkin(7) + (-t64 * t67 * t69 + qJDD(3)) * pkin(3) + t91;
t81 = -t67 * g(3) + t64 * t39;
t16 = -t62 * t69 * pkin(3) + t50 * pkin(7) - qJD(3) * t54 + t81;
t93 = t63 * t15 + t66 * t16;
t27 = -t43 * mrSges(6,1) + t44 * mrSges(6,2);
t86 = t43 * t27 + t23 * mrSges(6,3) + m(6) * (-t42 * pkin(4) + t23 * qJ(5) + 0.2e1 * qJD(5) * t43 - t60 * t35 + t93);
t59 = qJDD(3) + qJDD(4);
t79 = t66 * t15 - t63 * t16;
t85 = t60 * t33 + t59 * mrSges(6,1) + m(6) * (-0.2e1 * qJD(5) * t44 + (t43 * t60 - t24) * qJ(5) + (t43 * t44 + t59) * pkin(4) + t79);
t49 = (mrSges(4,1) * t64 + mrSges(4,2) * t67) * qJD(1);
t28 = -t43 * mrSges(5,1) + t44 * mrSges(5,2);
t6 = m(5) * t79 + t59 * mrSges(5,1) + t60 * t34 + (-t28 - t27) * t44 + (-mrSges(5,3) - mrSges(6,3)) * t24 + t85;
t7 = m(5) * t93 + t23 * mrSges(5,3) + t43 * t28 + (-t37 - t36) * t60 + (-mrSges(5,2) - mrSges(6,2)) * t59 + t86;
t3 = m(4) * t91 + qJDD(3) * mrSges(4,1) - t51 * mrSges(4,3) + qJD(3) * t52 - t49 * t89 + t66 * t6 + t63 * t7;
t4 = m(4) * t81 - qJDD(3) * mrSges(4,2) + t50 * mrSges(4,3) - qJD(3) * t53 - t49 * t90 - t63 * t6 + t66 * t7;
t83 = -t64 * t3 + t67 * t4;
t74 = -m(3) * (-qJDD(1) * pkin(1) + t73) - t67 * t3 - t64 * t4;
t5 = m(2) * t78 + t94 * qJDD(1) - t96 * t69 - t101;
t1 = m(2) * t82 + t96 * qJDD(1) + t94 * t69 + t74;
t2 = [-m(1) * g(1) - t65 * t1 + t68 * t5, t5, -m(3) * g(3) + t83, t4, t7, -t59 * mrSges(6,2) - t60 * t36 + t86; -m(1) * g(2) + t68 * t1 + t65 * t5, t1, -t69 * mrSges(3,2) - qJDD(1) * mrSges(3,3) + t101, t3, t6, -t24 * mrSges(6,3) - t44 * t27 + t85; (-m(1) + t100) * g(3) + t83, t100 * g(3) + t83, qJDD(1) * mrSges(3,2) - t69 * mrSges(3,3) - t74, t70, t72, -t23 * mrSges(6,1) - t43 * t33 + t84;];
f_new = t2;
