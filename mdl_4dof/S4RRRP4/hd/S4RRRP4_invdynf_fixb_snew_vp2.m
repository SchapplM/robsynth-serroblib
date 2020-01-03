% Calculate vector of cutting forces with Newton-Euler
% S4RRRP4
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x5]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP4_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP4_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP4_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:11
% EndTime: 2019-12-31 17:15:12
% DurationCPUTime: 0.50s
% Computational Cost: add. (3117->121), mult. (6588->156), div. (0->0), fcn. (3917->6), ass. (0->56)
t55 = sin(qJ(2));
t58 = cos(qJ(2));
t73 = qJD(1) * qJD(2);
t44 = t55 * qJDD(1) + t58 * t73;
t45 = t58 * qJDD(1) - t55 * t73;
t75 = qJD(1) * t55;
t46 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t75;
t74 = qJD(1) * t58;
t47 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t74;
t60 = qJD(1) ^ 2;
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t39 = (t54 * t58 + t55 * t57) * qJD(1);
t23 = -t39 * qJD(3) - t54 * t44 + t57 * t45;
t38 = (-t54 * t55 + t57 * t58) * qJD(1);
t24 = t38 * qJD(3) + t57 * t44 + t54 * t45;
t52 = qJD(2) + qJD(3);
t31 = -t52 * mrSges(5,2) + t38 * mrSges(5,3);
t32 = -t52 * mrSges(4,2) + t38 * mrSges(4,3);
t35 = t52 * mrSges(4,1) - t39 * mrSges(4,3);
t48 = qJD(2) * pkin(2) - pkin(6) * t75;
t53 = t58 ^ 2;
t56 = sin(qJ(1));
t59 = cos(qJ(1));
t69 = t56 * g(1) - t59 * g(2);
t63 = -qJDD(1) * pkin(1) - t69;
t61 = -t45 * pkin(2) + t48 * t75 + (-pkin(6) * t53 - pkin(5)) * t60 + t63;
t33 = t52 * pkin(3) - t39 * qJ(4);
t34 = t52 * mrSges(5,1) - t39 * mrSges(5,3);
t37 = t38 ^ 2;
t70 = m(5) * (-t23 * pkin(3) - t37 * qJ(4) + t39 * t33 + qJDD(4) + t61) + t24 * mrSges(5,2) + t39 * t34;
t62 = m(4) * t61 + t24 * mrSges(4,2) - (mrSges(4,1) + mrSges(5,1)) * t23 + t39 * t35 - (t32 + t31) * t38 + t70;
t86 = t62 + (t55 * t46 - t58 * t47) * qJD(1) - t45 * mrSges(3,1) + t44 * mrSges(3,2) + m(3) * (-t60 * pkin(5) + t63);
t43 = (-mrSges(3,1) * t58 + mrSges(3,2) * t55) * qJD(1);
t28 = -t38 * mrSges(5,1) + t39 * mrSges(5,2);
t29 = -t38 * mrSges(4,1) + t39 * mrSges(4,2);
t51 = qJDD(2) + qJDD(3);
t66 = -t59 * g(1) - t56 * g(2);
t41 = -t60 * pkin(1) + qJDD(1) * pkin(5) + t66;
t79 = t55 * t41;
t82 = pkin(2) * t60;
t16 = qJDD(2) * pkin(2) - t44 * pkin(6) - t79 + (pkin(6) * t73 + t55 * t82 - g(3)) * t58;
t68 = -t55 * g(3) + t58 * t41;
t17 = t45 * pkin(6) - qJD(2) * t48 - t53 * t82 + t68;
t67 = t57 * t16 - t54 * t17;
t72 = t52 * t31 + t51 * mrSges(5,1) + m(5) * (-0.2e1 * qJD(4) * t39 + (t38 * t52 - t24) * qJ(4) + (t38 * t39 + t51) * pkin(3) + t67);
t7 = m(4) * t67 + t51 * mrSges(4,1) + t52 * t32 + (-t29 - t28) * t39 + (-mrSges(4,3) - mrSges(5,3)) * t24 + t72;
t77 = t54 * t16 + t57 * t17;
t71 = m(5) * (-t37 * pkin(3) + t23 * qJ(4) + 0.2e1 * qJD(4) * t38 - t52 * t33 + t77) + t38 * t28 + t23 * mrSges(5,3);
t8 = m(4) * t77 + t23 * mrSges(4,3) + t38 * t29 + (-t35 - t34) * t52 + (-mrSges(4,2) - mrSges(5,2)) * t51 + t71;
t4 = m(3) * (-t58 * g(3) - t79) - t44 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t43 * t75 + qJD(2) * t47 + t54 * t8 + t57 * t7;
t5 = m(3) * t68 - qJDD(2) * mrSges(3,2) + t45 * mrSges(3,3) - qJD(2) * t46 + t43 * t74 - t54 * t7 + t57 * t8;
t84 = t58 * t4 + t55 * t5;
t6 = m(2) * t69 + qJDD(1) * mrSges(2,1) - t60 * mrSges(2,2) - t86;
t1 = m(2) * t66 - t60 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t55 * t4 + t58 * t5;
t2 = [-m(1) * g(1) + t59 * t1 - t56 * t6, t1, t5, t8, -t51 * mrSges(5,2) - t52 * t34 + t71; -m(1) * g(2) + t56 * t1 + t59 * t6, t6, t4, t7, -t24 * mrSges(5,3) - t39 * t28 + t72; (-m(1) - m(2)) * g(3) + t84, -m(2) * g(3) + t84, t86, t62, -t23 * mrSges(5,1) - t38 * t31 + t70;];
f_new = t2;
