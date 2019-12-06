% Calculate vector of cutting forces with Newton-Euler
% S5PRRPP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRPP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:47
% EndTime: 2019-12-05 16:08:49
% DurationCPUTime: 0.68s
% Computational Cost: add. (5349->127), mult. (11436->167), div. (0->0), fcn. (6941->8), ass. (0->64)
t64 = sin(qJ(3));
t66 = cos(qJ(3));
t84 = qJD(2) * qJD(3);
t80 = t66 * t84;
t48 = t64 * qJDD(2) + t80;
t49 = t66 * qJDD(2) - t64 * t84;
t86 = qJD(2) * t64;
t53 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t86;
t85 = qJD(2) * t66;
t54 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t85;
t69 = qJD(2) ^ 2;
t62 = sin(pkin(8));
t87 = cos(pkin(8));
t28 = t62 * t48 - t87 * t49;
t29 = t87 * t48 + t62 * t49;
t39 = t62 * t86 - t87 * t85;
t34 = -qJD(3) * mrSges(5,2) - t39 * mrSges(5,3);
t40 = (t62 * t66 + t87 * t64) * qJD(2);
t35 = qJD(3) * mrSges(5,1) - t40 * mrSges(5,3);
t52 = qJD(3) * pkin(3) - qJ(4) * t86;
t60 = t66 ^ 2;
t63 = sin(pkin(7));
t88 = cos(pkin(7));
t51 = -t88 * g(1) - t63 * g(2);
t61 = -g(3) + qJDD(1);
t65 = sin(qJ(2));
t67 = cos(qJ(2));
t78 = -t65 * t51 + t67 * t61;
t75 = -qJDD(2) * pkin(2) - t78;
t71 = -t49 * pkin(3) + qJDD(4) + t52 * t86 + (-qJ(4) * t60 - pkin(6)) * t69 + t75;
t36 = -qJD(3) * mrSges(6,1) + t40 * mrSges(6,2);
t37 = -t39 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t73 = t29 * mrSges(6,3) + t40 * t36 - m(6) * (-0.2e1 * qJD(5) * t40 + (qJD(3) * t39 - t29) * qJ(5) + (qJD(3) * t40 + t28) * pkin(4) + t71) - t39 * t37 - t28 * mrSges(6,1);
t72 = m(5) * t71 + t28 * mrSges(5,1) + t29 * mrSges(5,2) + t39 * t34 + t40 * t35 - t73;
t95 = (t64 * t53 - t66 * t54) * qJD(2) + m(4) * (-t69 * pkin(6) + t75) - t49 * mrSges(4,1) + t48 * mrSges(4,2) + t72;
t94 = -2 * qJD(4);
t23 = t39 * pkin(4) - t40 * qJ(5);
t68 = qJD(3) ^ 2;
t89 = t67 * t51 + t65 * t61;
t32 = -t69 * pkin(2) + qJDD(2) * pkin(6) + t89;
t50 = t63 * g(1) - t88 * g(2);
t79 = -t64 * t32 - t66 * t50;
t17 = (-t48 + t80) * qJ(4) + (t64 * t66 * t69 + qJDD(3)) * pkin(3) + t79;
t90 = t66 * t32 - t64 * t50;
t18 = -t60 * t69 * pkin(3) + t49 * qJ(4) - qJD(3) * t52 + t90;
t74 = t87 * t17 - t62 * t18;
t93 = m(6) * (-qJDD(3) * pkin(4) - t68 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t23) * t40 - t74);
t92 = -mrSges(5,3) - mrSges(6,2);
t24 = t39 * mrSges(6,1) - t40 * mrSges(6,3);
t91 = -t39 * mrSges(5,1) - t40 * mrSges(5,2) - t24;
t10 = m(5) * t74 - t93 + (m(5) * t94 + t91) * t40 + t92 * t29 + (mrSges(5,1) + mrSges(6,1)) * qJDD(3) + (t34 + t37) * qJD(3);
t47 = (-mrSges(4,1) * t66 + mrSges(4,2) * t64) * qJD(2);
t81 = t62 * t17 + t87 * t18 + t39 * t94;
t82 = m(6) * (-t68 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t39 * t23 + t81) + qJD(3) * t36 + qJDD(3) * mrSges(6,3);
t9 = m(5) * t81 - qJDD(3) * mrSges(5,2) - qJD(3) * t35 + t92 * t28 + t91 * t39 + t82;
t5 = m(4) * t79 + qJDD(3) * mrSges(4,1) - t48 * mrSges(4,3) + qJD(3) * t54 + t87 * t10 - t47 * t86 + t62 * t9;
t6 = m(4) * t90 - qJDD(3) * mrSges(4,2) + t49 * mrSges(4,3) - qJD(3) * t53 - t62 * t10 + t47 * t85 + t87 * t9;
t3 = m(3) * t89 - t69 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t64 * t5 + t66 * t6;
t8 = m(3) * t78 + qJDD(2) * mrSges(3,1) - t69 * mrSges(3,2) - t95;
t83 = m(2) * t61 + t65 * t3 + t67 * t8;
t77 = -t66 * t5 - t64 * t6;
t4 = (m(2) + m(3)) * t50 + t77;
t1 = m(2) * t51 + t67 * t3 - t65 * t8;
t2 = [-m(1) * g(1) + t88 * t1 - t63 * t4, t1, t3, t6, t9, -t28 * mrSges(6,2) - t39 * t24 + t82; -m(1) * g(2) + t63 * t1 + t88 * t4, t4, t8, t5, t10, -t73; -m(1) * g(3) + t83, t83, -m(3) * t50 - t77, t95, t72, -qJDD(3) * mrSges(6,1) + t29 * mrSges(6,2) - qJD(3) * t37 + t40 * t24 + t93;];
f_new = t2;
