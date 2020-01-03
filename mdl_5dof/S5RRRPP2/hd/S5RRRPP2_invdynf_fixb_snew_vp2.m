% Calculate vector of cutting forces with Newton-Euler
% S5RRRPP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:31
% EndTime: 2019-12-31 20:51:32
% DurationCPUTime: 0.53s
% Computational Cost: add. (3646->140), mult. (4602->165), div. (0->0), fcn. (1994->6), ass. (0->64)
t57 = qJDD(1) + qJDD(2);
t63 = sin(qJ(3));
t58 = qJD(1) + qJD(2);
t66 = cos(qJ(3));
t86 = qJD(3) * t66;
t83 = t58 * t86;
t35 = t63 * t57 + t83;
t97 = t58 * t63;
t36 = -qJD(3) * t97 + t66 * t57;
t47 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t97;
t48 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t97;
t65 = sin(qJ(1));
t68 = cos(qJ(1));
t81 = t65 * g(1) - t68 * g(2);
t43 = qJDD(1) * pkin(1) + t81;
t70 = qJD(1) ^ 2;
t76 = -t68 * g(1) - t65 * g(2);
t44 = -t70 * pkin(1) + t76;
t64 = sin(qJ(2));
t67 = cos(qJ(2));
t90 = t67 * t43 - t64 * t44;
t78 = t57 * pkin(2) + t90;
t74 = -t35 * qJ(4) - t78;
t101 = 2 * qJD(4);
t45 = -qJD(3) * pkin(4) - qJ(5) * t97;
t46 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t97;
t96 = t58 * t66;
t49 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t96;
t56 = t58 ^ 2;
t62 = t66 ^ 2;
t77 = m(6) * (qJDD(5) + (-qJ(5) * t62 + pkin(7)) * t56 + (pkin(3) + pkin(4)) * t36 + (qJ(4) * t86 + (-pkin(3) * qJD(3) + t101 + t45) * t63) * t58 - t74) + t35 * mrSges(6,2) + t36 * mrSges(6,1) + t46 * t97 + t49 * t96;
t98 = t56 * pkin(7);
t75 = m(5) * (-t36 * pkin(3) - t98 + (-0.2e1 * qJD(4) * t63 + (pkin(3) * t63 - qJ(4) * t66) * qJD(3)) * t58 + t74) - t36 * mrSges(5,1) - t77;
t51 = mrSges(5,2) * t96 + qJD(3) * mrSges(5,3);
t87 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t96 + t51;
t104 = (-t87 * t66 + (t47 - t48) * t63) * t58 + (mrSges(4,2) - mrSges(5,3)) * t35 + m(4) * (-t78 - t98) - t36 * mrSges(4,1) + t75;
t100 = -m(2) - m(3);
t32 = (-mrSges(5,1) * t66 - mrSges(5,3) * t63) * t58;
t31 = (-pkin(3) * t66 - qJ(4) * t63) * t58;
t69 = qJD(3) ^ 2;
t89 = t64 * t43 + t67 * t44;
t21 = -t56 * pkin(2) + t57 * pkin(7) + t89;
t80 = -t63 * g(3) + t66 * t21;
t72 = -t69 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t101 + t31 * t96 + t80;
t84 = -0.2e1 * qJD(5) * t58;
t79 = m(6) * (-t62 * t56 * pkin(4) - t36 * qJ(5) + qJD(3) * t45 + t66 * t84 + t72) + qJD(3) * t46 + qJDD(3) * mrSges(6,2) - t36 * mrSges(6,3);
t73 = m(5) * t72 + qJDD(3) * mrSges(5,3) + qJD(3) * t48 + t32 * t96 + t79;
t33 = (mrSges(6,1) * t66 + mrSges(6,2) * t63) * t58;
t91 = -t33 + (-mrSges(4,1) * t66 + mrSges(4,2) * t63) * t58;
t93 = mrSges(4,3) + mrSges(5,2);
t7 = m(4) * t80 - qJDD(3) * mrSges(4,2) - qJD(3) * t47 + t93 * t36 + t91 * t96 + t73;
t92 = -t66 * g(3) - t63 * t21;
t17 = -qJDD(3) * pkin(3) - t69 * qJ(4) + t31 * t97 + qJDD(4) - t92;
t12 = m(6) * (t63 * t84 + (-t35 + t83) * qJ(5) + (-t56 * t63 * t66 - qJDD(3)) * pkin(4) + t17);
t82 = m(5) * t17 + t12;
t95 = -mrSges(5,1) - mrSges(6,1);
t8 = m(4) * t92 + (-t32 - t91) * t97 + (mrSges(6,3) - t93) * t35 + (mrSges(4,1) - t95) * qJDD(3) + (t49 + t87) * qJD(3) - t82;
t99 = t63 * t7 + t66 * t8;
t85 = t33 * t96;
t4 = m(3) * t90 + t57 * mrSges(3,1) - t56 * mrSges(3,2) - t104;
t3 = m(3) * t89 - t56 * mrSges(3,1) - t57 * mrSges(3,2) - t63 * t8 + t66 * t7;
t2 = m(2) * t76 - t70 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t67 * t3 - t64 * t4;
t1 = m(2) * t81 + qJDD(1) * mrSges(2,1) - t70 * mrSges(2,2) + t64 * t3 + t67 * t4;
t5 = [-m(1) * g(1) - t65 * t1 + t68 * t2, t2, t3, t7, t36 * mrSges(5,2) + t73 - t85, t79 - t85; -m(1) * g(2) + t68 * t1 + t65 * t2, t1, t4, t8, -t35 * mrSges(5,3) + (-t63 * t48 - t66 * t51) * t58 + t75, -qJDD(3) * mrSges(6,1) - t35 * mrSges(6,3) - qJD(3) * t49 - t33 * t97 + t12; (-m(1) + t100) * g(3) + t99, t100 * g(3) + t99, -m(3) * g(3) + t99, t104, (t32 - t33) * t97 + (mrSges(5,2) - mrSges(6,3)) * t35 + t95 * qJDD(3) + (-t49 - t51) * qJD(3) + t82, t77;];
f_new = t5;
