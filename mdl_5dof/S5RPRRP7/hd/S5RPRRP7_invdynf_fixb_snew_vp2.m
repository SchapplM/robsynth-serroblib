% Calculate vector of cutting forces with Newton-Euler
% S5RPRRP7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRP7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP7_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP7_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP7_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:36
% EndTime: 2019-12-31 18:44:38
% DurationCPUTime: 0.66s
% Computational Cost: add. (5970->133), mult. (11202->167), div. (0->0), fcn. (6351->8), ass. (0->67)
t62 = sin(qJ(4));
t63 = sin(qJ(3));
t84 = qJD(1) * t63;
t92 = cos(qJ(4));
t42 = -t92 * qJD(3) + t62 * t84;
t65 = cos(qJ(3));
t82 = qJD(1) * qJD(3);
t76 = t65 * t82;
t48 = t63 * qJDD(1) + t76;
t27 = -t42 * qJD(4) + t62 * qJDD(3) + t92 * t48;
t83 = t65 * qJD(1);
t53 = qJD(4) - t83;
t32 = -t53 * mrSges(5,2) - t42 * mrSges(5,3);
t35 = -t42 * mrSges(6,2) + t53 * mrSges(6,3);
t77 = t63 * t82;
t49 = t65 * qJDD(1) - t77;
t41 = qJDD(4) - t49;
t43 = t62 * qJD(3) + t92 * t84;
t68 = qJD(1) ^ 2;
t64 = sin(qJ(1));
t66 = cos(qJ(1));
t78 = t64 * g(1) - t66 * g(2);
t44 = qJDD(1) * pkin(1) + t78;
t73 = -t66 * g(1) - t64 * g(2);
t46 = -t68 * pkin(1) + t73;
t60 = sin(pkin(8));
t61 = cos(pkin(8));
t74 = t61 * t44 - t60 * t46;
t22 = -qJDD(1) * pkin(2) - t68 * pkin(6) - t74;
t16 = (-t48 - t76) * pkin(7) + (-t49 + t77) * pkin(3) + t22;
t47 = (-pkin(3) * t65 - pkin(7) * t63) * qJD(1);
t67 = qJD(3) ^ 2;
t85 = t60 * t44 + t61 * t46;
t23 = -t68 * pkin(2) + qJDD(1) * pkin(6) + t85;
t59 = -g(3) + qJDD(2);
t88 = t65 * t23 + t63 * t59;
t19 = -t67 * pkin(3) + qJDD(3) * pkin(7) + t47 * t83 + t88;
t71 = t92 * t16 - t62 * t19;
t30 = t42 * mrSges(6,1) - t43 * mrSges(6,3);
t87 = -t42 * mrSges(5,1) - t43 * mrSges(5,2) - t30;
t90 = -mrSges(5,3) - mrSges(6,2);
t29 = t42 * pkin(4) - t43 * qJ(5);
t52 = t53 ^ 2;
t93 = m(6) * (-t41 * pkin(4) - t52 * qJ(5) + t43 * t29 + qJDD(5) - t71);
t10 = m(5) * t71 - t93 + (t32 + t35) * t53 + t87 * t43 + (mrSges(5,1) + mrSges(6,1)) * t41 + t90 * t27;
t50 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t84;
t51 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t83;
t26 = t43 * qJD(4) - t92 * qJDD(3) + t62 * t48;
t33 = t53 * mrSges(5,1) - t43 * mrSges(5,3);
t34 = -t53 * mrSges(6,1) + t43 * mrSges(6,2);
t89 = t62 * t16 + t92 * t19;
t80 = m(6) * (-t52 * pkin(4) + t41 * qJ(5) + 0.2e1 * qJD(5) * t53 - t42 * t29 + t89) + t53 * t34 + t41 * mrSges(6,3);
t9 = m(5) * t89 - t41 * mrSges(5,2) + t90 * t26 - t53 * t33 + t87 * t42 + t80;
t95 = m(4) * t22 - t49 * mrSges(4,1) + t48 * mrSges(4,2) + (t50 * t63 - t51 * t65) * qJD(1) + t92 * t10 + t62 * t9;
t75 = -t63 * t23 + t65 * t59;
t18 = -qJDD(3) * pkin(3) - t67 * pkin(7) + t47 * t84 - t75;
t79 = m(6) * (-0.2e1 * qJD(5) * t43 + (t42 * t53 - t27) * qJ(5) + (t43 * t53 + t26) * pkin(4) + t18) + t26 * mrSges(6,1) + t42 * t35;
t94 = m(5) * t18 + t26 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t27 + t42 * t32 + (t33 - t34) * t43 + t79;
t45 = (-mrSges(4,1) * t65 + mrSges(4,2) * t63) * qJD(1);
t6 = m(4) * t88 - qJDD(3) * mrSges(4,2) + t49 * mrSges(4,3) - qJD(3) * t50 - t62 * t10 + t45 * t83 + t92 * t9;
t8 = m(4) * t75 + qJDD(3) * mrSges(4,1) - t48 * mrSges(4,3) + qJD(3) * t51 - t45 * t84 - t94;
t81 = m(3) * t59 + t63 * t6 + t65 * t8;
t4 = m(3) * t74 + qJDD(1) * mrSges(3,1) - t68 * mrSges(3,2) - t95;
t3 = m(3) * t85 - t68 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t65 * t6 - t63 * t8;
t2 = m(2) * t73 - t68 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t61 * t3 - t60 * t4;
t1 = m(2) * t78 + qJDD(1) * mrSges(2,1) - t68 * mrSges(2,2) + t60 * t3 + t61 * t4;
t5 = [-m(1) * g(1) - t64 * t1 + t66 * t2, t2, t3, t6, t9, -t26 * mrSges(6,2) - t42 * t30 + t80; -m(1) * g(2) + t66 * t1 + t64 * t2, t1, t4, t8, t10, -t27 * mrSges(6,3) - t43 * t34 + t79; (-m(1) - m(2)) * g(3) + t81, -m(2) * g(3) + t81, t81, t95, t94, -t41 * mrSges(6,1) + t27 * mrSges(6,2) + t43 * t30 - t53 * t35 + t93;];
f_new = t5;
