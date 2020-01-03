% Calculate vector of cutting forces with Newton-Euler
% S4RPRP6
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RPRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP6_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_invdynf_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP6_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP6_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:00
% EndTime: 2019-12-31 16:46:00
% DurationCPUTime: 0.27s
% Computational Cost: add. (847->93), mult. (1614->112), div. (0->0), fcn. (625->4), ass. (0->43)
t44 = qJD(1) ^ 2;
t40 = sin(qJ(3));
t42 = cos(qJ(3));
t61 = qJD(1) * qJD(3);
t26 = -t40 * qJDD(1) - t42 * t61;
t54 = t40 * t61;
t27 = t42 * qJDD(1) - t54;
t63 = qJD(1) * t40;
t30 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t63;
t62 = qJD(1) * t42;
t33 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t62;
t41 = sin(qJ(1));
t43 = cos(qJ(1));
t51 = -t43 * g(1) - t41 * g(2);
t49 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t51;
t29 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t63;
t31 = qJD(3) * pkin(3) - qJ(4) * t62;
t32 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t62;
t39 = t40 ^ 2;
t70 = -pkin(1) - pkin(5);
t53 = t29 * t63 + t32 * t62 + t27 * mrSges(5,2) + m(5) * (t31 * t62 - t26 * pkin(3) + qJDD(4) + (-qJ(4) * t39 + t70) * t44 + t49);
t45 = -(mrSges(4,1) + mrSges(5,1)) * t26 + m(4) * (t70 * t44 + t49) + t30 * t63 + t33 * t62 + t27 * mrSges(4,2) + t53;
t72 = m(3) * (t44 * pkin(1) - t49) - t45;
t71 = -m(2) - m(3);
t56 = t41 * g(1) - t43 * g(2);
t46 = -t44 * qJ(2) + qJDD(2) - t56;
t14 = t70 * qJDD(1) + t46;
t55 = -t42 * g(3) + t40 * t14;
t58 = -0.2e1 * qJD(1) * qJD(4);
t68 = t26 * mrSges(5,3) + m(5) * (-t39 * t44 * pkin(3) + t26 * qJ(4) - qJD(3) * t31 + t40 * t58 + t55);
t67 = mrSges(2,1) - mrSges(3,2);
t65 = -mrSges(2,2) + mrSges(3,3);
t64 = t40 * g(3) + t42 * t14;
t59 = qJD(3) * t29 + qJDD(3) * mrSges(5,1) + m(5) * (t42 * t58 + (-t27 - t54) * qJ(4) + (-t40 * t42 * t44 + qJDD(3)) * pkin(3) + t64);
t24 = (mrSges(5,1) * t40 + mrSges(5,2) * t42) * qJD(1);
t52 = qJD(1) * (-t24 - (mrSges(4,1) * t40 + mrSges(4,2) * t42) * qJD(1));
t4 = m(4) * t64 + qJDD(3) * mrSges(4,1) + qJD(3) * t30 + (-mrSges(4,3) - mrSges(5,3)) * t27 + t42 * t52 + t59;
t5 = m(4) * t55 + t26 * mrSges(4,3) + (-mrSges(4,2) - mrSges(5,2)) * qJDD(3) + (-t33 - t32) * qJD(3) + t40 * t52 + t68;
t57 = -t40 * t4 + t42 * t5;
t48 = -m(3) * (-qJDD(1) * pkin(1) + t46) - t42 * t4 - t40 * t5;
t2 = m(2) * t51 + t65 * qJDD(1) - t67 * t44 - t72;
t1 = m(2) * t56 + t67 * qJDD(1) + t65 * t44 + t48;
t3 = [-m(1) * g(1) - t41 * t1 + t43 * t2, t2, -m(3) * g(3) + t57, t5, -qJDD(3) * mrSges(5,2) - qJD(3) * t32 - t24 * t63 + t68; -m(1) * g(2) + t43 * t1 + t41 * t2, t1, -t44 * mrSges(3,2) - qJDD(1) * mrSges(3,3) + t72, t4, -t27 * mrSges(5,3) - t24 * t62 + t59; (-m(1) + t71) * g(3) + t57, t71 * g(3) + t57, qJDD(1) * mrSges(3,2) - t44 * mrSges(3,3) - t48, t45, -t26 * mrSges(5,1) + t53;];
f_new = t3;
