% Calculate vector of cutting forces with Newton-Euler
% S5RPRRP13
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRP13_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP13_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP13_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP13_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:33
% EndTime: 2019-12-31 18:58:34
% DurationCPUTime: 0.51s
% Computational Cost: add. (3856->134), mult. (7245->160), div. (0->0), fcn. (3892->6), ass. (0->66)
t60 = sin(qJ(1));
t62 = cos(qJ(1));
t74 = -t62 * g(1) - t60 * g(2);
t98 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t74;
t59 = sin(qJ(3));
t61 = cos(qJ(3));
t46 = (pkin(3) * t59 - pkin(7) * t61) * qJD(1);
t63 = qJD(3) ^ 2;
t64 = qJD(1) ^ 2;
t78 = t60 * g(1) - t62 * g(2);
t69 = -t64 * qJ(2) + qJDD(2) - t78;
t95 = -pkin(1) - pkin(6);
t34 = t95 * qJDD(1) + t69;
t73 = t59 * g(3) + t61 * t34;
t85 = qJD(1) * t61;
t17 = -qJDD(3) * pkin(3) - t63 * pkin(7) + t46 * t85 - t73;
t58 = sin(qJ(4));
t93 = cos(qJ(4));
t44 = t58 * qJD(3) + t93 * t85;
t83 = qJD(1) * qJD(3);
t76 = t59 * t83;
t48 = t61 * qJDD(1) - t76;
t21 = t44 * qJD(4) - t93 * qJDD(3) + t58 * t48;
t43 = -t93 * qJD(3) + t58 * t85;
t22 = -t43 * qJD(4) + t58 * qJDD(3) + t93 * t48;
t84 = t59 * qJD(1);
t52 = qJD(4) + t84;
t27 = -t52 * mrSges(5,2) - t43 * mrSges(5,3);
t28 = t52 * mrSges(5,1) - t44 * mrSges(5,3);
t29 = -t52 * mrSges(6,1) + t44 * mrSges(6,2);
t30 = -t43 * mrSges(6,2) + t52 * mrSges(6,3);
t80 = m(6) * (-0.2e1 * qJD(5) * t44 + (t43 * t52 - t22) * qJ(5) + (t44 * t52 + t21) * pkin(4) + t17) + t21 * mrSges(6,1) + t43 * t30;
t97 = m(5) * t17 + t21 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t22 + t43 * t27 + (t28 - t29) * t44 + t80;
t96 = -m(2) - m(3);
t24 = t43 * pkin(4) - t44 * qJ(5);
t75 = t61 * t83;
t47 = -t59 * qJDD(1) - t75;
t42 = qJDD(4) - t47;
t51 = t52 ^ 2;
t65 = t95 * t64 - t98;
t15 = (-t48 + t76) * pkin(7) + (-t47 + t75) * pkin(3) + t65;
t77 = -t61 * g(3) + t59 * t34;
t18 = -t63 * pkin(3) + qJDD(3) * pkin(7) - t46 * t84 + t77;
t71 = t93 * t15 - t58 * t18;
t94 = m(6) * (-t42 * pkin(4) - t51 * qJ(5) + t44 * t24 + qJDD(5) - t71);
t92 = (mrSges(2,1) - mrSges(3,2));
t91 = -mrSges(2,2) + mrSges(3,3);
t89 = -mrSges(5,3) - mrSges(6,2);
t88 = t58 * t15 + t93 * t18;
t25 = t43 * mrSges(6,1) - t44 * mrSges(6,3);
t87 = -t43 * mrSges(5,1) - t44 * mrSges(5,2) - t25;
t81 = m(6) * (-t51 * pkin(4) + t42 * qJ(5) + 0.2e1 * qJD(5) * t52 - t43 * t24 + t88) + t52 * t29 + t42 * mrSges(6,3);
t45 = (mrSges(4,1) * t59 + mrSges(4,2) * t61) * qJD(1);
t50 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t85;
t7 = m(5) * t88 - t42 * mrSges(5,2) + t89 * t21 - t52 * t28 + t87 * t43 + t81;
t9 = m(5) * t71 - t94 + (t27 + t30) * t52 + t87 * t44 + (mrSges(5,1) + mrSges(6,1)) * t42 + t89 * t22;
t4 = m(4) * t77 - qJDD(3) * mrSges(4,2) + t47 * mrSges(4,3) - qJD(3) * t50 - t45 * t84 - t58 * t9 + t93 * t7;
t49 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t84;
t5 = m(4) * t73 + qJDD(3) * mrSges(4,1) - t48 * mrSges(4,3) + qJD(3) * t49 - t45 * t85 - t97;
t79 = t61 * t4 - t59 * t5;
t70 = -m(3) * (-qJDD(1) * pkin(1) + t69) - t59 * t4 - t61 * t5;
t68 = m(4) * t65 - t47 * mrSges(4,1) + t48 * mrSges(4,2) + t49 * t84 + t50 * t85 + t58 * t7 + t93 * t9;
t66 = -m(3) * (t64 * pkin(1) + t98) + t68;
t2 = m(2) * t74 + t91 * qJDD(1) - (t92 * t64) + t66;
t1 = m(2) * t78 + t92 * qJDD(1) + t91 * t64 + t70;
t3 = [-m(1) * g(1) - t60 * t1 + t62 * t2, t2, -m(3) * g(3) + t79, t4, t7, -t21 * mrSges(6,2) - t43 * t25 + t81; -m(1) * g(2) + t62 * t1 + t60 * t2, t1, -(t64 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t66, t5, t9, -t22 * mrSges(6,3) - t44 * t29 + t80; (-m(1) + t96) * g(3) + t79, t96 * g(3) + t79, qJDD(1) * mrSges(3,2) - t64 * mrSges(3,3) - t70, t68, t97, -t42 * mrSges(6,1) + t22 * mrSges(6,2) + t44 * t25 - t52 * t30 + t94;];
f_new = t3;
