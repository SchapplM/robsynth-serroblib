% Calculate vector of cutting forces with Newton-Euler
% S5RPRRR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRR12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR12_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR12_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR12_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:31
% EndTime: 2019-12-31 19:12:34
% DurationCPUTime: 0.81s
% Computational Cost: add. (8224->138), mult. (16112->179), div. (0->0), fcn. (10048->8), ass. (0->73)
t68 = sin(qJ(1));
t72 = cos(qJ(1));
t85 = -t72 * g(1) - t68 * g(2);
t81 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t85;
t66 = sin(qJ(4));
t67 = sin(qJ(3));
t70 = cos(qJ(4));
t71 = cos(qJ(3));
t45 = (t66 * t71 + t67 * t70) * qJD(1);
t99 = -m(2) - m(3);
t98 = -pkin(1) - pkin(6);
t97 = (mrSges(2,1) - mrSges(3,2));
t96 = -mrSges(2,2) + mrSges(3,3);
t91 = qJD(1) * qJD(3);
t86 = t67 * t91;
t53 = t71 * qJDD(1) - t86;
t73 = qJD(1) ^ 2;
t88 = t68 * g(1) - t72 * g(2);
t79 = -t73 * qJ(2) + qJDD(2) - t88;
t41 = t98 * qJDD(1) + t79;
t94 = t67 * g(3) + t71 * t41;
t19 = (-t53 - t86) * pkin(7) + (-t67 * t71 * t73 + qJDD(3)) * pkin(3) + t94;
t52 = -t67 * qJDD(1) - t71 * t91;
t92 = qJD(1) * t71;
t56 = qJD(3) * pkin(3) - pkin(7) * t92;
t64 = t67 ^ 2;
t87 = -t71 * g(3) + t67 * t41;
t20 = -t64 * t73 * pkin(3) + t52 * pkin(7) - qJD(3) * t56 + t87;
t95 = t66 * t19 + t70 * t20;
t93 = qJD(1) * t67;
t51 = (mrSges(4,1) * t67 + mrSges(4,2) * t71) * qJD(1);
t54 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t93;
t46 = (-t66 * t67 + t70 * t71) * qJD(1);
t26 = -t46 * qJD(4) + t70 * t52 - t66 * t53;
t27 = -t45 * qJD(4) + t66 * t52 + t70 * t53;
t62 = qJD(3) + qJD(4);
t75 = -t52 * pkin(3) + t56 * t92 + (-pkin(7) * t64 + t98) * t73 + t81;
t12 = (t45 * t62 - t27) * pkin(8) + (t46 * t62 - t26) * pkin(4) + t75;
t32 = t45 * pkin(4) - t46 * pkin(8);
t60 = t62 ^ 2;
t61 = qJDD(3) + qJDD(4);
t14 = -t60 * pkin(4) + t61 * pkin(8) - t45 * t32 + t95;
t65 = sin(qJ(5));
t69 = cos(qJ(5));
t33 = -t65 * t46 + t69 * t62;
t16 = t33 * qJD(5) + t69 * t27 + t65 * t61;
t34 = t69 * t46 + t65 * t62;
t21 = -t33 * mrSges(6,1) + t34 * mrSges(6,2);
t25 = qJDD(5) - t26;
t44 = qJD(5) + t45;
t28 = -t44 * mrSges(6,2) + t33 * mrSges(6,3);
t10 = m(6) * (t69 * t12 - t65 * t14) - t16 * mrSges(6,3) + t25 * mrSges(6,1) - t34 * t21 + t44 * t28;
t15 = -t34 * qJD(5) - t65 * t27 + t69 * t61;
t29 = t44 * mrSges(6,1) - t34 * mrSges(6,3);
t11 = m(6) * (t65 * t12 + t69 * t14) + t15 * mrSges(6,3) - t25 * mrSges(6,2) + t33 * t21 - t44 * t29;
t31 = t45 * mrSges(5,1) + t46 * mrSges(5,2);
t39 = t62 * mrSges(5,1) - t46 * mrSges(5,3);
t6 = m(5) * t95 - t61 * mrSges(5,2) + t26 * mrSges(5,3) - t65 * t10 + t69 * t11 - t45 * t31 - t62 * t39;
t38 = -t62 * mrSges(5,2) - t45 * mrSges(5,3);
t84 = t70 * t19 - t66 * t20;
t77 = m(6) * (-t61 * pkin(4) - t60 * pkin(8) + t46 * t32 - t84) - t15 * mrSges(6,1) + t16 * mrSges(6,2) - t33 * t28 + t34 * t29;
t7 = m(5) * t84 + t61 * mrSges(5,1) - t27 * mrSges(5,3) - t46 * t31 + t62 * t38 - t77;
t3 = m(4) * t94 + qJDD(3) * mrSges(4,1) - t53 * mrSges(4,3) + qJD(3) * t54 - t51 * t92 + t66 * t6 + t70 * t7;
t55 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t92;
t4 = m(4) * t87 - qJDD(3) * mrSges(4,2) + t52 * mrSges(4,3) - qJD(3) * t55 - t51 * t93 + t70 * t6 - t66 * t7;
t89 = -t67 * t3 + t71 * t4;
t80 = -m(3) * (-qJDD(1) * pkin(1) + t79) - t71 * t3 - t67 * t4;
t78 = m(5) * t75 - t26 * mrSges(5,1) + t27 * mrSges(5,2) + t69 * t10 + t65 * t11 + t45 * t38 + t46 * t39;
t76 = -t52 * mrSges(4,1) + m(4) * (t98 * t73 + t81) + t54 * t93 + t55 * t92 + t53 * mrSges(4,2) + t78;
t74 = -m(3) * (t73 * pkin(1) - t81) + t76;
t5 = m(2) * t85 + t96 * qJDD(1) - (t97 * t73) + t74;
t1 = m(2) * t88 + t97 * qJDD(1) + t96 * t73 + t80;
t2 = [-m(1) * g(1) - t68 * t1 + t72 * t5, t5, -m(3) * g(3) + t89, t4, t6, t11; -m(1) * g(2) + t72 * t1 + t68 * t5, t1, -(t73 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t74, t3, t7, t10; (-m(1) + t99) * g(3) + t89, t99 * g(3) + t89, qJDD(1) * mrSges(3,2) - t73 * mrSges(3,3) - t80, t76, t78, t77;];
f_new = t2;
