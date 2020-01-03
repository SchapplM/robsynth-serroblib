% Calculate vector of cutting forces with Newton-Euler
% S5RPRRP12
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRP12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP12_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP12_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP12_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:17
% EndTime: 2019-12-31 18:56:19
% DurationCPUTime: 0.52s
% Computational Cost: add. (3911->134), mult. (7425->160), div. (0->0), fcn. (4037->6), ass. (0->65)
t61 = sin(qJ(1));
t64 = cos(qJ(1));
t75 = -t64 * g(1) - t61 * g(2);
t97 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t75;
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t50 = (pkin(3) * t60 - pkin(7) * t63) * qJD(1);
t65 = qJD(3) ^ 2;
t66 = qJD(1) ^ 2;
t80 = t61 * g(1) - t64 * g(2);
t71 = -t66 * qJ(2) + qJDD(2) - t80;
t94 = -pkin(1) - pkin(6);
t37 = t94 * qJDD(1) + t71;
t74 = t60 * g(3) + t63 * t37;
t88 = qJD(1) * t63;
t18 = -qJDD(3) * pkin(3) - t65 * pkin(7) + t50 * t88 - t74;
t59 = sin(qJ(4));
t62 = cos(qJ(4));
t48 = t59 * qJD(3) + t62 * t88;
t86 = qJD(1) * qJD(3);
t78 = t60 * t86;
t52 = t63 * qJDD(1) - t78;
t24 = -t48 * qJD(4) + t62 * qJDD(3) - t59 * t52;
t47 = t62 * qJD(3) - t59 * t88;
t25 = t47 * qJD(4) + t59 * qJDD(3) + t62 * t52;
t87 = t60 * qJD(1);
t55 = qJD(4) + t87;
t29 = -t55 * mrSges(6,2) + t47 * mrSges(6,3);
t30 = -t55 * mrSges(5,2) + t47 * mrSges(5,3);
t33 = t55 * mrSges(5,1) - t48 * mrSges(5,3);
t31 = t55 * pkin(4) - t48 * qJ(5);
t32 = t55 * mrSges(6,1) - t48 * mrSges(6,3);
t45 = t47 ^ 2;
t82 = m(6) * (-t24 * pkin(4) - t45 * qJ(5) + t48 * t31 + qJDD(5) + t18) + t25 * mrSges(6,2) + t48 * t32;
t96 = m(5) * t18 + t25 * mrSges(5,2) - (mrSges(5,1) + mrSges(6,1)) * t24 + t48 * t33 - (t30 + t29) * t47 + t82;
t95 = -m(2) - m(3);
t93 = (mrSges(2,1) - mrSges(3,2));
t91 = -mrSges(2,2) + mrSges(3,3);
t77 = t63 * t86;
t51 = -t60 * qJDD(1) - t77;
t67 = t94 * t66 - t97;
t16 = (-t52 + t78) * pkin(7) + (-t51 + t77) * pkin(3) + t67;
t79 = -t63 * g(3) + t60 * t37;
t19 = -t65 * pkin(3) + qJDD(3) * pkin(7) - t50 * t87 + t79;
t90 = t59 * t16 + t62 * t19;
t46 = qJDD(4) - t51;
t76 = t62 * t16 - t59 * t19;
t84 = m(6) * (-0.2e1 * qJD(5) * t48 + (t47 * t55 - t25) * qJ(5) + (t47 * t48 + t46) * pkin(4) + t76) + t55 * t29 + t46 * mrSges(6,1);
t27 = -t47 * mrSges(6,1) + t48 * mrSges(6,2);
t83 = m(6) * (-t45 * pkin(4) + t24 * qJ(5) + 0.2e1 * qJD(5) * t47 - t55 * t31 + t90) + t47 * t27 + t24 * mrSges(6,3);
t49 = (mrSges(4,1) * t60 + mrSges(4,2) * t63) * qJD(1);
t54 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t88;
t28 = -t47 * mrSges(5,1) + t48 * mrSges(5,2);
t6 = m(5) * t76 + t46 * mrSges(5,1) + t55 * t30 + (-t28 - t27) * t48 + (-mrSges(5,3) - mrSges(6,3)) * t25 + t84;
t8 = m(5) * t90 + t24 * mrSges(5,3) + t47 * t28 + (-t33 - t32) * t55 + (-mrSges(5,2) - mrSges(6,2)) * t46 + t83;
t4 = m(4) * t79 - qJDD(3) * mrSges(4,2) + t51 * mrSges(4,3) - qJD(3) * t54 - t49 * t87 - t59 * t6 + t62 * t8;
t53 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t87;
t9 = m(4) * t74 + qJDD(3) * mrSges(4,1) - t52 * mrSges(4,3) + qJD(3) * t53 - t49 * t88 - t96;
t81 = t63 * t4 - t60 * t9;
t72 = -m(3) * (-qJDD(1) * pkin(1) + t71) - t60 * t4 - t63 * t9;
t70 = m(4) * t67 - t51 * mrSges(4,1) + t52 * mrSges(4,2) + t53 * t87 + t54 * t88 + t59 * t8 + t62 * t6;
t68 = -m(3) * (t66 * pkin(1) + t97) + t70;
t2 = m(2) * t75 + t91 * qJDD(1) - (t93 * t66) + t68;
t1 = m(2) * t80 + t93 * qJDD(1) + t91 * t66 + t72;
t3 = [-m(1) * g(1) - t61 * t1 + t64 * t2, t2, -m(3) * g(3) + t81, t4, t8, -t46 * mrSges(6,2) - t55 * t32 + t83; -m(1) * g(2) + t64 * t1 + t61 * t2, t1, -(t66 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t68, t9, t6, -t25 * mrSges(6,3) - t48 * t27 + t84; (-m(1) + t95) * g(3) + t81, t95 * g(3) + t81, qJDD(1) * mrSges(3,2) - t66 * mrSges(3,3) - t72, t70, t96, -t24 * mrSges(6,1) - t47 * t29 + t82;];
f_new = t3;
