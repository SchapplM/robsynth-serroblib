% Calculate vector of cutting forces with Newton-Euler
% S5RPPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRR12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR12_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR12_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR12_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:06:59
% EndTime: 2019-12-31 18:07:01
% DurationCPUTime: 0.74s
% Computational Cost: add. (6385->125), mult. (14160->158), div. (0->0), fcn. (9156->8), ass. (0->72)
t68 = qJD(1) ^ 2;
t63 = sin(qJ(1));
t66 = cos(qJ(1));
t85 = t63 * g(1) - t66 * g(2);
t75 = -t68 * qJ(2) + qJDD(2) - t85;
t95 = -pkin(1) - qJ(3);
t102 = -(2 * qJD(1) * qJD(3)) + t95 * qJDD(1) + t75;
t59 = sin(pkin(8));
t56 = t59 ^ 2;
t60 = cos(pkin(8));
t93 = t60 ^ 2 + t56;
t84 = t93 * mrSges(4,3);
t82 = -t66 * g(1) - t63 * g(2);
t101 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t82;
t62 = sin(qJ(4));
t65 = cos(qJ(4));
t80 = t59 * t65 + t60 * t62;
t45 = t80 * qJD(1);
t79 = -t59 * t62 + t60 * t65;
t46 = t79 * qJD(1);
t91 = t46 * qJD(4);
t31 = -t80 * qJDD(1) - t91;
t100 = -m(2) - m(3);
t99 = pkin(3) * t68;
t98 = mrSges(4,2) * t60;
t97 = mrSges(2,1) - mrSges(3,2);
t96 = -mrSges(2,2) + mrSges(3,3);
t88 = t59 * g(3) + t102 * t60;
t20 = (-pkin(6) * qJDD(1) - t59 * t99) * t60 + t88;
t83 = -t60 * g(3) + t102 * t59;
t90 = qJDD(1) * t59;
t21 = -pkin(6) * t90 - t56 * t99 + t83;
t94 = t62 * t20 + t65 * t21;
t92 = t45 * qJD(4);
t30 = t45 * pkin(4) - t46 * pkin(7);
t67 = qJD(4) ^ 2;
t13 = -t67 * pkin(4) + qJDD(4) * pkin(7) - t45 * t30 + t94;
t32 = t79 * qJDD(1) - t92;
t74 = qJDD(3) + t101;
t69 = pkin(3) * t90 + (-t93 * pkin(6) + t95) * t68 + t74;
t14 = (-t32 + t92) * pkin(7) + (-t31 + t91) * pkin(4) + t69;
t61 = sin(qJ(5));
t64 = cos(qJ(5));
t34 = t64 * qJD(4) - t61 * t46;
t16 = t34 * qJD(5) + t61 * qJDD(4) + t64 * t32;
t35 = t61 * qJD(4) + t64 * t46;
t18 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t43 = qJD(5) + t45;
t22 = -t43 * mrSges(6,2) + t34 * mrSges(6,3);
t29 = qJDD(5) - t31;
t10 = m(6) * (-t61 * t13 + t64 * t14) - t16 * mrSges(6,3) + t29 * mrSges(6,1) - t35 * t18 + t43 * t22;
t15 = -t35 * qJD(5) + t64 * qJDD(4) - t61 * t32;
t23 = t43 * mrSges(6,1) - t35 * mrSges(6,3);
t11 = m(6) * (t64 * t13 + t61 * t14) + t15 * mrSges(6,3) - t29 * mrSges(6,2) + t34 * t18 - t43 * t23;
t27 = t45 * mrSges(5,1) + t46 * mrSges(5,2);
t40 = qJD(4) * mrSges(5,1) - t46 * mrSges(5,3);
t6 = m(5) * t94 - qJDD(4) * mrSges(5,2) + t31 * mrSges(5,3) - qJD(4) * t40 - t61 * t10 + t64 * t11 - t45 * t27;
t39 = -qJD(4) * mrSges(5,2) - t45 * mrSges(5,3);
t81 = t65 * t20 - t62 * t21;
t71 = m(6) * (-qJDD(4) * pkin(4) - t67 * pkin(7) + t46 * t30 - t81) - t15 * mrSges(6,1) + t16 * mrSges(6,2) - t34 * t22 + t35 * t23;
t7 = m(5) * t81 + qJDD(4) * mrSges(5,1) - t32 * mrSges(5,3) + qJD(4) * t39 - t46 * t27 - t71;
t77 = -qJDD(1) * mrSges(4,3) - t68 * (mrSges(4,1) * t59 + t98);
t3 = m(4) * t88 + t62 * t6 + t77 * t60 + t65 * t7;
t4 = m(4) * t83 + t77 * t59 + t65 * t6 - t62 * t7;
t86 = -t59 * t3 + t60 * t4;
t76 = -m(3) * (-qJDD(1) * pkin(1) + t75) - t60 * t3 - t59 * t4;
t73 = -m(5) * t69 + t31 * mrSges(5,1) - t32 * mrSges(5,2) - t64 * t10 - t61 * t11 - t45 * t39 - t46 * t40;
t72 = m(4) * (t95 * t68 + t74) + mrSges(4,1) * t90 + qJDD(1) * t98 - t73;
t70 = m(3) * (t68 * pkin(1) - t101) - t72;
t5 = -t70 + m(2) * t82 + t96 * qJDD(1) + (-t84 - t97) * t68;
t1 = m(2) * t85 + t97 * qJDD(1) + t96 * t68 + t76;
t2 = [-m(1) * g(1) - t63 * t1 + t66 * t5, t5, -m(3) * g(3) + t86, t4, t6, t11; -m(1) * g(2) + t66 * t1 + t63 * t5, t1, -qJDD(1) * mrSges(3,3) + t70 + (-mrSges(3,2) + t84) * t68, t3, t7, t10; (-m(1) + t100) * g(3) + t86, t100 * g(3) + t86, qJDD(1) * mrSges(3,2) - t68 * mrSges(3,3) - t76, -t68 * t84 + t72, -t73, t71;];
f_new = t2;
