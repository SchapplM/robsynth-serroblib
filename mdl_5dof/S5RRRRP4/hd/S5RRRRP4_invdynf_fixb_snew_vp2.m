% Calculate vector of cutting forces with Newton-Euler
% S5RRRRP4
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:43
% EndTime: 2019-12-31 21:50:45
% DurationCPUTime: 0.74s
% Computational Cost: add. (10008->136), mult. (12678->173), div. (0->0), fcn. (7367->8), ass. (0->69)
t59 = qJDD(1) + qJDD(2);
t64 = sin(qJ(3));
t67 = cos(qJ(3));
t61 = qJD(1) + qJD(2);
t83 = qJD(3) * t61;
t43 = t64 * t59 + t67 * t83;
t44 = t67 * t59 - t64 * t83;
t90 = t61 * t64;
t50 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t90;
t89 = t61 * t67;
t51 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t89;
t57 = t61 ^ 2;
t63 = sin(qJ(4));
t91 = cos(qJ(4));
t39 = (t63 * t67 + t91 * t64) * t61;
t23 = t39 * qJD(4) + t63 * t43 - t91 * t44;
t38 = t63 * t90 - t91 * t89;
t24 = -t38 * qJD(4) + t91 * t43 + t63 * t44;
t60 = qJD(3) + qJD(4);
t34 = -t60 * mrSges(5,2) - t38 * mrSges(5,3);
t35 = t60 * mrSges(5,1) - t39 * mrSges(5,3);
t52 = qJD(3) * pkin(3) - pkin(8) * t90;
t62 = t67 ^ 2;
t66 = sin(qJ(1));
t69 = cos(qJ(1));
t81 = t66 * g(1) - t69 * g(2);
t48 = qJDD(1) * pkin(1) + t81;
t70 = qJD(1) ^ 2;
t78 = -t69 * g(1) - t66 * g(2);
t49 = -t70 * pkin(1) + t78;
t65 = sin(qJ(2));
t68 = cos(qJ(2));
t79 = t68 * t48 - t65 * t49;
t76 = -t59 * pkin(2) - t79;
t73 = -t44 * pkin(3) + t52 * t90 + (-pkin(8) * t62 - pkin(7)) * t57 + t76;
t36 = -t60 * mrSges(6,1) + t39 * mrSges(6,2);
t37 = -t38 * mrSges(6,2) + t60 * mrSges(6,3);
t74 = t24 * mrSges(6,3) + t39 * t36 - m(6) * (-0.2e1 * qJD(5) * t39 + (t38 * t60 - t24) * qJ(5) + (t39 * t60 + t23) * pkin(4) + t73) - t23 * mrSges(6,1) - t38 * t37;
t72 = m(5) * t73 + t23 * mrSges(5,1) + t24 * mrSges(5,2) + t38 * t34 + t39 * t35 - t74;
t96 = (t64 * t50 - t67 * t51) * t61 + m(4) * (-t57 * pkin(7) + t76) - t44 * mrSges(4,1) + t43 * mrSges(4,2) + t72;
t95 = -m(2) - m(3);
t58 = qJDD(3) + qJDD(4);
t84 = t65 * t48 + t68 * t49;
t32 = -t57 * pkin(2) + t59 * pkin(7) + t84;
t88 = t64 * t32;
t92 = pkin(3) * t57;
t17 = qJDD(3) * pkin(3) - t43 * pkin(8) - t88 + (pkin(8) * t83 + t64 * t92 - g(3)) * t67;
t80 = -t64 * g(3) + t67 * t32;
t18 = t44 * pkin(8) - qJD(3) * t52 - t62 * t92 + t80;
t75 = t91 * t17 - t63 * t18;
t29 = t38 * mrSges(6,1) - t39 * mrSges(6,3);
t85 = -t38 * mrSges(5,1) - t39 * mrSges(5,2) - t29;
t87 = -mrSges(5,3) - mrSges(6,2);
t28 = t38 * pkin(4) - t39 * qJ(5);
t56 = t60 ^ 2;
t93 = m(6) * (-t58 * pkin(4) - t56 * qJ(5) + t39 * t28 + qJDD(5) - t75);
t10 = m(5) * t75 - t93 + (t34 + t37) * t60 + (mrSges(5,1) + mrSges(6,1)) * t58 + t85 * t39 + t87 * t24;
t42 = (-mrSges(4,1) * t67 + mrSges(4,2) * t64) * t61;
t86 = t63 * t17 + t91 * t18;
t82 = m(6) * (-t56 * pkin(4) + t58 * qJ(5) + 0.2e1 * qJD(5) * t60 - t38 * t28 + t86) + t60 * t36 + t58 * mrSges(6,3);
t9 = m(5) * t86 - t58 * mrSges(5,2) + t87 * t23 - t60 * t35 + t85 * t38 + t82;
t6 = m(4) * (-t67 * g(3) - t88) - t43 * mrSges(4,3) + qJDD(3) * mrSges(4,1) - t42 * t90 + qJD(3) * t51 + t63 * t9 + t91 * t10;
t7 = m(4) * t80 - qJDD(3) * mrSges(4,2) + t44 * mrSges(4,3) - qJD(3) * t50 - t63 * t10 + t42 * t89 + t91 * t9;
t94 = t67 * t6 + t64 * t7;
t8 = m(3) * t79 + t59 * mrSges(3,1) - t57 * mrSges(3,2) - t96;
t3 = m(3) * t84 - t57 * mrSges(3,1) - t59 * mrSges(3,2) - t64 * t6 + t67 * t7;
t2 = m(2) * t78 - t70 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t68 * t3 - t65 * t8;
t1 = m(2) * t81 + qJDD(1) * mrSges(2,1) - t70 * mrSges(2,2) + t65 * t3 + t68 * t8;
t4 = [-m(1) * g(1) - t66 * t1 + t69 * t2, t2, t3, t7, t9, -t23 * mrSges(6,2) - t38 * t29 + t82; -m(1) * g(2) + t69 * t1 + t66 * t2, t1, t8, t6, t10, -t74; (-m(1) + t95) * g(3) + t94, t95 * g(3) + t94, -m(3) * g(3) + t94, t96, t72, -t58 * mrSges(6,1) + t24 * mrSges(6,2) + t39 * t29 - t60 * t37 + t93;];
f_new = t4;
