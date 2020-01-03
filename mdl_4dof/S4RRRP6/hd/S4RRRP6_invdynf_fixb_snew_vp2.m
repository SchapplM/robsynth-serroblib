% Calculate vector of cutting forces with Newton-Euler
% S4RRRP6
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP6_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP6_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP6_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:10
% EndTime: 2019-12-31 17:18:12
% DurationCPUTime: 0.42s
% Computational Cost: add. (2807->119), mult. (5546->150), div. (0->0), fcn. (3187->6), ass. (0->57)
t59 = qJD(1) ^ 2;
t54 = sin(qJ(1));
t57 = cos(qJ(1));
t68 = t54 * g(1) - t57 * g(2);
t36 = -qJDD(1) * pkin(1) - t59 * pkin(5) - t68;
t53 = sin(qJ(2));
t56 = cos(qJ(2));
t72 = qJD(1) * qJD(2);
t65 = t56 * t72;
t45 = t53 * qJDD(1) + t65;
t66 = t53 * t72;
t46 = t56 * qJDD(1) - t66;
t74 = qJD(1) * t53;
t47 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t74;
t73 = t56 * qJD(1);
t48 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t73;
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t41 = t55 * qJD(2) - t52 * t74;
t24 = t41 * qJD(3) + t52 * qJDD(2) + t55 * t45;
t42 = t52 * qJD(2) + t55 * t74;
t26 = -t41 * mrSges(5,1) + t42 * mrSges(5,2);
t27 = -t41 * mrSges(4,1) + t42 * mrSges(4,2);
t49 = qJD(3) - t73;
t29 = -t49 * mrSges(4,2) + t41 * mrSges(4,3);
t40 = qJDD(3) - t46;
t15 = (-t45 - t65) * pkin(6) + (-t46 + t66) * pkin(2) + t36;
t44 = (-pkin(2) * t56 - pkin(6) * t53) * qJD(1);
t58 = qJD(2) ^ 2;
t63 = -t57 * g(1) - t54 * g(2);
t37 = -t59 * pkin(1) + qJDD(1) * pkin(5) + t63;
t67 = -t53 * g(3) + t56 * t37;
t18 = -t58 * pkin(2) + qJDD(2) * pkin(6) + t44 * t73 + t67;
t64 = t55 * t15 - t52 * t18;
t28 = -t49 * mrSges(5,2) + t41 * mrSges(5,3);
t71 = t49 * t28 + t40 * mrSges(5,1) + m(5) * (-0.2e1 * qJD(4) * t42 + (t41 * t49 - t24) * qJ(4) + (t41 * t42 + t40) * pkin(3) + t64);
t5 = m(4) * t64 + t40 * mrSges(4,1) + t49 * t29 + (-t27 - t26) * t42 + (-mrSges(4,3) - mrSges(5,3)) * t24 + t71;
t23 = -t42 * qJD(3) + t55 * qJDD(2) - t52 * t45;
t31 = t49 * mrSges(5,1) - t42 * mrSges(5,3);
t32 = t49 * mrSges(4,1) - t42 * mrSges(4,3);
t30 = t49 * pkin(3) - t42 * qJ(4);
t39 = t41 ^ 2;
t77 = t52 * t15 + t55 * t18;
t70 = m(5) * (-t39 * pkin(3) + t23 * qJ(4) + 0.2e1 * qJD(4) * t41 - t49 * t30 + t77) + t41 * t26 + t23 * mrSges(5,3);
t6 = m(4) * t77 + t23 * mrSges(4,3) + t41 * t27 + (-t32 - t31) * t49 + (-mrSges(4,2) - mrSges(5,2)) * t40 + t70;
t81 = m(3) * t36 - t46 * mrSges(3,1) + t45 * mrSges(3,2) + t55 * t5 + t52 * t6 + (t53 * t47 - t56 * t48) * qJD(1);
t75 = -t56 * g(3) - t53 * t37;
t17 = -qJDD(2) * pkin(2) - t58 * pkin(6) + t44 * t74 - t75;
t69 = m(5) * (-t23 * pkin(3) - t39 * qJ(4) + t42 * t30 + qJDD(4) + t17) + t24 * mrSges(5,2) + t42 * t31;
t80 = m(4) * t17 + t24 * mrSges(4,2) - (mrSges(4,1) + mrSges(5,1)) * t23 + t42 * t32 - (t29 + t28) * t41 + t69;
t43 = (-mrSges(3,1) * t56 + mrSges(3,2) * t53) * qJD(1);
t4 = m(3) * t67 - qJDD(2) * mrSges(3,2) + t46 * mrSges(3,3) - qJD(2) * t47 + t43 * t73 - t52 * t5 + t55 * t6;
t8 = m(3) * t75 + qJDD(2) * mrSges(3,1) - t45 * mrSges(3,3) + qJD(2) * t48 - t43 * t74 - t80;
t79 = t53 * t4 + t56 * t8;
t2 = m(2) * t68 + qJDD(1) * mrSges(2,1) - t59 * mrSges(2,2) - t81;
t1 = m(2) * t63 - t59 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t56 * t4 - t53 * t8;
t3 = [-m(1) * g(1) + t57 * t1 - t54 * t2, t1, t4, t6, -t40 * mrSges(5,2) - t49 * t31 + t70; -m(1) * g(2) + t54 * t1 + t57 * t2, t2, t8, t5, -t24 * mrSges(5,3) - t42 * t26 + t71; (-m(1) - m(2)) * g(3) + t79, -m(2) * g(3) + t79, t81, t80, -t23 * mrSges(5,1) - t41 * t28 + t69;];
f_new = t3;
