% Calculate vector of cutting forces with Newton-Euler
% S5PRPRP2
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:34
% EndTime: 2019-12-05 15:30:36
% DurationCPUTime: 0.55s
% Computational Cost: add. (3855->114), mult. (8293->155), div. (0->0), fcn. (4988->8), ass. (0->69)
t59 = sin(pkin(8));
t55 = t59 ^ 2;
t61 = cos(pkin(8));
t100 = (t61 ^ 2 + t55) * mrSges(4,3);
t73 = -pkin(3) * t61 - pkin(6) * t59;
t44 = t73 * qJD(2);
t67 = qJD(2) ^ 2;
t60 = sin(pkin(7));
t62 = cos(pkin(7));
t45 = t60 * g(1) - t62 * g(2);
t46 = -t62 * g(1) - t60 * g(2);
t64 = sin(qJ(2));
t66 = cos(qJ(2));
t94 = t64 * t45 + t66 * t46;
t24 = -t67 * pkin(2) + qJDD(2) * qJ(3) + t94;
t58 = -g(3) + qJDD(1);
t78 = -t59 * t24 + t61 * t58;
t81 = 0.2e1 * qJD(2) * qJD(3);
t92 = qJD(2) * t59;
t16 = t44 * t92 + t59 * t81 - t78;
t63 = sin(qJ(4));
t65 = cos(qJ(4));
t89 = qJD(2) * qJD(4);
t37 = (-qJDD(2) * t63 - t65 * t89) * t59;
t38 = (qJDD(2) * t65 - t63 * t89) * t59;
t91 = t61 * qJD(2);
t49 = qJD(4) - t91;
t84 = t63 * t92;
t30 = -t49 * mrSges(6,2) - mrSges(6,3) * t84;
t80 = qJ(5) * t92;
t32 = t49 * pkin(4) - t65 * t80;
t83 = t65 * t92;
t33 = t49 * mrSges(6,1) - mrSges(6,3) * t83;
t98 = t55 * t67;
t87 = t63 ^ 2 * t98;
t76 = m(6) * (-t37 * pkin(4) - qJ(5) * t87 + t32 * t83 + qJDD(5) + t16) + t30 * t84 + t33 * t83 + t38 * mrSges(6,2);
t99 = -m(5) * t16 - t38 * mrSges(5,2) + (mrSges(5,1) + mrSges(6,1)) * t37 - t76;
t75 = -0.2e1 * qJD(5) * t92;
t85 = t59 * t58 + (t24 + t81) * t61;
t17 = t44 * t91 + t85;
t77 = t66 * t45 - t64 * t46;
t69 = -t67 * qJ(3) + qJDD(3) - t77;
t20 = (-pkin(2) + t73) * qJDD(2) + t69;
t95 = t65 * t17 + t63 * t20;
t96 = m(6) * (-pkin(4) * t87 + t37 * qJ(5) - t49 * t32 + t63 * t75 + t95) + t37 * mrSges(6,3);
t90 = qJDD(2) * mrSges(4,3);
t72 = -t61 * mrSges(4,1) + t59 * mrSges(4,2);
t43 = t72 * qJD(2);
t31 = -t49 * mrSges(5,2) - mrSges(5,3) * t84;
t34 = t49 * mrSges(5,1) - mrSges(5,3) * t83;
t71 = t31 * t63 + t34 * t65;
t10 = m(4) * t78 + (-t90 + (-0.2e1 * m(4) * qJD(3) - t43 - t71) * qJD(2)) * t59 + t99;
t19 = t65 * t20;
t48 = -t61 * qJDD(2) + qJDD(4);
t35 = (mrSges(6,1) * t63 + mrSges(6,2) * t65) * t92;
t74 = (-t35 - (mrSges(5,1) * t63 + mrSges(5,2) * t65) * t92) * t92;
t86 = m(6) * (t65 * t75 + t48 * pkin(4) - t38 * qJ(5) + t19 + (-pkin(4) * t65 * t98 - t49 * t80 - t17) * t63) + t49 * t30 + t48 * mrSges(6,1);
t7 = m(5) * (-t63 * t17 + t19) + t48 * mrSges(5,1) + t49 * t31 + (-mrSges(5,3) - mrSges(6,3)) * t38 + t65 * t74 + t86;
t8 = m(5) * t95 + t37 * mrSges(5,3) + (-t34 - t33) * t49 + (-mrSges(5,2) - mrSges(6,2)) * t48 + t63 * t74 + t96;
t6 = m(4) * t85 + t65 * t8 - t63 * t7 + (qJD(2) * t43 + t90) * t61;
t88 = m(3) * t58 + t61 * t10 + t59 * t6;
t82 = t35 * t92;
t79 = m(2) * t58 + t88;
t70 = m(4) * (-qJDD(2) * pkin(2) + t69) + t63 * t8 + t65 * t7;
t4 = m(3) * t77 + (-mrSges(3,2) + t100) * t67 + (mrSges(3,1) - t72) * qJDD(2) - t70;
t3 = m(3) * t94 - t67 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t59 * t10 + t61 * t6;
t2 = m(2) * t46 + t66 * t3 - t64 * t4;
t1 = m(2) * t45 + t64 * t3 + t66 * t4;
t5 = [-m(1) * g(1) - t60 * t1 + t62 * t2, t2, t3, t6, t8, -t48 * mrSges(6,2) - t49 * t33 - t63 * t82 + t96; -m(1) * g(2) + t62 * t1 + t60 * t2, t1, t4, t10, t7, -t38 * mrSges(6,3) - t65 * t82 + t86; -m(1) * g(3) + t79, t79, t88, t72 * qJDD(2) - t67 * t100 + t70, t71 * t92 - t99, -t37 * mrSges(6,1) + t76;];
f_new = t5;
