% Calculate vector of cutting forces with Newton-Euler
% S5RPRPP5
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
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:08
% EndTime: 2019-12-31 18:16:10
% DurationCPUTime: 0.44s
% Computational Cost: add. (1521->138), mult. (2934->158), div. (0->0), fcn. (1116->4), ass. (0->61)
t100 = -2 * qJD(1);
t65 = qJD(1) ^ 2;
t61 = sin(qJ(1));
t63 = cos(qJ(1));
t89 = -t63 * g(1) - t61 * g(2);
t19 = (t65 * pkin(1)) - qJDD(1) * qJ(2) + (qJD(2) * t100) - t89;
t99 = 2 * qJD(4);
t98 = -m(2) - m(3);
t97 = -pkin(3) - pkin(4);
t96 = pkin(4) * t65;
t95 = (mrSges(2,1) - mrSges(3,2));
t94 = -mrSges(5,1) - mrSges(6,1);
t93 = -mrSges(2,2) + mrSges(3,3);
t92 = -mrSges(4,3) - mrSges(5,2);
t60 = sin(qJ(3));
t62 = cos(qJ(3));
t34 = (mrSges(5,1) * t60 - mrSges(5,3) * t62) * qJD(1);
t35 = (-mrSges(6,1) * t60 + mrSges(6,2) * t62) * qJD(1);
t91 = t34 - t35;
t86 = qJD(1) * t60;
t41 = qJD(3) * mrSges(6,2) + mrSges(6,3) * t86;
t47 = -mrSges(5,2) * t86 + qJD(3) * mrSges(5,3);
t90 = -t41 - t47;
t88 = qJ(4) * t60;
t84 = qJD(1) * qJD(3);
t38 = t62 * qJDD(1) - t60 * t84;
t87 = t38 * qJ(4);
t85 = qJD(1) * t62;
t33 = (pkin(3) * t60 - qJ(4) * t62) * qJD(1);
t64 = qJD(3) ^ 2;
t74 = -t64 * qJ(4) + t33 * t85 + qJDD(4);
t80 = t61 * g(1) - t63 * g(2);
t70 = -t65 * qJ(2) + qJDD(2) - t80;
t18 = (-pkin(1) - pkin(6)) * qJDD(1) + t70;
t75 = t60 * g(3) + t62 * t18;
t9 = m(6) * (-t38 * qJ(5) + ((qJD(5) * t100) - t18) * t62 + t97 * qJDD(3) + (-qJ(5) * t84 + t62 * t96 - g(3)) * t60 + t74);
t82 = m(5) * (-qJDD(3) * pkin(3) + t74 - t75) + t9;
t36 = (mrSges(4,1) * t60 + mrSges(4,2) * t62) * qJD(1);
t37 = t60 * qJDD(1) + t62 * t84;
t45 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t85;
t46 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t85;
t79 = -t62 * g(3) + t60 * t18;
t69 = -t64 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t99 + t79;
t43 = -qJD(3) * pkin(4) - qJ(5) * t85;
t44 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t85;
t58 = t60 ^ 2;
t77 = t35 * t86 + t37 * mrSges(6,3) + qJD(3) * t44 + m(6) * (-t58 * t96 + t37 * qJ(5) + qJD(3) * t43 + ((2 * qJD(5)) - t33) * t86 + t69) + qJDD(3) * mrSges(6,2);
t71 = m(5) * (-t33 * t86 + t69) + qJD(3) * t46 + qJDD(3) * mrSges(5,3) + t77;
t4 = m(4) * t79 - qJDD(3) * mrSges(4,2) - qJD(3) * t45 + t92 * t37 + (-t34 - t36) * t86 + t71;
t42 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t86;
t5 = m(4) * t75 + (mrSges(6,3) + t92) * t38 + (mrSges(4,1) - t94) * qJDD(3) + (t42 - t90) * qJD(3) + (-t36 - t91) * t85 - t82;
t81 = t62 * t4 - t60 * t5;
t76 = m(6) * (t87 + qJDD(5) + (-qJ(5) * t58 + pkin(6)) * t65 + t97 * t37 + (-qJD(3) * t88 + (-pkin(3) * qJD(3) + t43 + t99) * t62) * qJD(1) + t19) + t38 * mrSges(6,2) - t37 * mrSges(6,1) + t44 * t85 - t41 * t86;
t73 = -m(3) * (-qJDD(1) * pkin(1) + t70) - t60 * t4 - t62 * t5;
t72 = -(t65 * pkin(6)) - t19;
t68 = -t38 * mrSges(5,3) - t46 * t85 + m(5) * (t37 * pkin(3) - t87 + (-0.2e1 * qJD(4) * t62 + (pkin(3) * t62 + t88) * qJD(3)) * qJD(1) + t72) + t47 * t86 + t37 * mrSges(5,1) - t76;
t67 = m(4) * t72 + t37 * mrSges(4,1) + t38 * mrSges(4,2) + t42 * t86 + t45 * t85 + t68;
t66 = -m(3) * t19 + t67;
t2 = m(2) * t89 + t93 * qJDD(1) - (t95 * t65) + t66;
t1 = m(2) * t80 + t95 * qJDD(1) + t93 * t65 + t73;
t3 = [-m(1) * g(1) - t61 * t1 + t63 * t2, t2, -m(3) * g(3) + t81, t4, -t37 * mrSges(5,2) - t34 * t86 + t71, t77; -m(1) * g(2) + t63 * t1 + t61 * t2, t1, -(t65 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t66, t5, t68, -qJDD(3) * mrSges(6,1) - t38 * mrSges(6,3) - qJD(3) * t41 - t35 * t85 + t9; (-m(1) + t98) * g(3) + t81, t98 * g(3) + t81, qJDD(1) * mrSges(3,2) - t65 * mrSges(3,3) - t73, t67, (mrSges(5,2) - mrSges(6,3)) * t38 + t94 * qJDD(3) + t90 * qJD(3) + t91 * t85 + t82, t76;];
f_new = t3;
