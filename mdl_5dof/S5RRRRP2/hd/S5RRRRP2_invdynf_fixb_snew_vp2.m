% Calculate vector of cutting forces with Newton-Euler
% S5RRRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:08
% EndTime: 2020-01-03 12:11:10
% DurationCPUTime: 0.80s
% Computational Cost: add. (10464->137), mult. (13296->174), div. (0->0), fcn. (7843->8), ass. (0->68)
t60 = qJDD(1) + qJDD(2);
t65 = sin(qJ(3));
t69 = cos(qJ(3));
t62 = qJD(1) + qJD(2);
t86 = qJD(3) * t62;
t46 = t65 * t60 + t69 * t86;
t47 = t69 * t60 - t65 * t86;
t93 = t62 * t65;
t53 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t93;
t92 = t62 * t69;
t54 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t92;
t58 = t62 ^ 2;
t64 = sin(qJ(4));
t68 = cos(qJ(4));
t43 = (t64 * t69 + t65 * t68) * t62;
t26 = -t43 * qJD(4) - t64 * t46 + t68 * t47;
t42 = (-t64 * t65 + t68 * t69) * t62;
t27 = t42 * qJD(4) + t68 * t46 + t64 * t47;
t61 = qJD(3) + qJD(4);
t36 = -t61 * mrSges(6,2) + t42 * mrSges(6,3);
t37 = -t61 * mrSges(5,2) + t42 * mrSges(5,3);
t40 = t61 * mrSges(5,1) - t43 * mrSges(5,3);
t55 = qJD(3) * pkin(3) - pkin(8) * t93;
t63 = t69 ^ 2;
t67 = sin(qJ(1));
t71 = cos(qJ(1));
t78 = -t71 * g(2) - t67 * g(3);
t51 = qJDD(1) * pkin(1) + t78;
t72 = qJD(1) ^ 2;
t81 = -t67 * g(2) + t71 * g(3);
t52 = -t72 * pkin(1) + t81;
t66 = sin(qJ(2));
t70 = cos(qJ(2));
t79 = t70 * t51 - t66 * t52;
t75 = -t60 * pkin(2) - t79;
t73 = -t47 * pkin(3) + t55 * t93 + (-pkin(8) * t63 - pkin(7)) * t58 + t75;
t38 = t61 * pkin(4) - t43 * qJ(5);
t39 = t61 * mrSges(6,1) - t43 * mrSges(6,3);
t41 = t42 ^ 2;
t83 = m(6) * (-t26 * pkin(4) - t41 * qJ(5) + t43 * t38 + qJDD(5) + t73) + t27 * mrSges(6,2) + t43 * t39;
t74 = m(5) * t73 + t27 * mrSges(5,2) - (mrSges(5,1) + mrSges(6,1)) * t26 + t43 * t40 - (t37 + t36) * t42 + t83;
t101 = (t65 * t53 - t69 * t54) * t62 - t47 * mrSges(4,1) + t46 * mrSges(4,2) + m(4) * (-t58 * pkin(7) + t75) + t74;
t99 = -m(2) - m(3);
t32 = -t42 * mrSges(5,1) + t43 * mrSges(5,2);
t59 = qJDD(3) + qJDD(4);
t31 = -t42 * mrSges(6,1) + t43 * mrSges(6,2);
t87 = t66 * t51 + t70 * t52;
t34 = -t58 * pkin(2) + t60 * pkin(7) + t87;
t91 = t65 * t34;
t96 = pkin(3) * t58;
t18 = qJDD(3) * pkin(3) - t46 * pkin(8) - t91 + (pkin(8) * t86 + t65 * t96 - g(1)) * t69;
t82 = -t65 * g(1) + t69 * t34;
t19 = t47 * pkin(8) - qJD(3) * t55 - t63 * t96 + t82;
t89 = t64 * t18 + t68 * t19;
t84 = m(6) * (-t41 * pkin(4) + t26 * qJ(5) + 0.2e1 * qJD(5) * t42 - t61 * t38 + t89) + t42 * t31 + t26 * mrSges(6,3);
t10 = m(5) * t89 + t26 * mrSges(5,3) + t42 * t32 + (-t40 - t39) * t61 + (-mrSges(5,2) - mrSges(6,2)) * t59 + t84;
t45 = (-mrSges(4,1) * t69 + mrSges(4,2) * t65) * t62;
t80 = t68 * t18 - t64 * t19;
t85 = m(6) * (-0.2e1 * qJD(5) * t43 + (t42 * t61 - t27) * qJ(5) + (t42 * t43 + t59) * pkin(4) + t80) + t61 * t36 + t59 * mrSges(6,1);
t9 = m(5) * t80 + t59 * mrSges(5,1) + t61 * t37 + (-t32 - t31) * t43 + (-mrSges(5,3) - mrSges(6,3)) * t27 + t85;
t6 = m(4) * (-t69 * g(1) - t91) - t46 * mrSges(4,3) + qJDD(3) * mrSges(4,1) - t45 * t93 + qJD(3) * t54 + t64 * t10 + t68 * t9;
t7 = m(4) * t82 - qJDD(3) * mrSges(4,2) + t47 * mrSges(4,3) - qJD(3) * t53 + t68 * t10 + t45 * t92 - t64 * t9;
t98 = t69 * t6 + t65 * t7;
t8 = m(3) * t79 + t60 * mrSges(3,1) - t58 * mrSges(3,2) - t101;
t3 = m(3) * t87 - t58 * mrSges(3,1) - t60 * mrSges(3,2) - t65 * t6 + t69 * t7;
t2 = m(2) * t78 + qJDD(1) * mrSges(2,1) - t72 * mrSges(2,2) + t66 * t3 + t70 * t8;
t1 = m(2) * t81 - t72 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t70 * t3 - t66 * t8;
t4 = [(-m(1) + t99) * g(1) + t98, t1, t3, t7, t10, -t59 * mrSges(6,2) - t61 * t39 + t84; -m(1) * g(2) + t67 * t1 + t71 * t2, t2, t8, t6, t9, -t27 * mrSges(6,3) - t43 * t31 + t85; -m(1) * g(3) - t71 * t1 + t67 * t2, t99 * g(1) + t98, -m(3) * g(1) + t98, t101, t74, -t26 * mrSges(6,1) - t42 * t36 + t83;];
f_new = t4;
