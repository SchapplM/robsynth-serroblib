% Calculate vector of cutting forces with Newton-Euler
% S5RPPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:23
% EndTime: 2019-12-31 17:46:24
% DurationCPUTime: 0.47s
% Computational Cost: add. (4412->101), mult. (8377->128), div. (0->0), fcn. (3995->8), ass. (0->58)
t54 = qJD(1) ^ 2;
t48 = cos(pkin(8));
t41 = t48 ^ 2;
t46 = sin(pkin(8));
t73 = t46 ^ 2 + t41;
t80 = t73 * mrSges(5,3);
t79 = -m(2) - m(3);
t78 = -pkin(1) - pkin(2);
t77 = pkin(4) * t54;
t76 = mrSges(2,1) + mrSges(3,1);
t51 = sin(qJ(1));
t53 = cos(qJ(1));
t67 = -t53 * g(1) - t51 * g(2);
t60 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t67;
t25 = t78 * t54 + t60;
t69 = t51 * g(1) - t53 * g(2);
t57 = -t54 * qJ(2) + qJDD(2) - t69;
t29 = t78 * qJDD(1) + t57;
t47 = sin(pkin(7));
t49 = cos(pkin(7));
t75 = t49 * t25 + t47 * t29;
t44 = g(3) + qJDD(3);
t71 = qJD(1) * qJD(4);
t74 = t48 * t44 + 0.2e1 * t46 * t71;
t72 = qJDD(1) * t48;
t16 = -t54 * pkin(3) - qJDD(1) * qJ(4) + t75;
t70 = t46 * t44 + (t16 - (2 * t71)) * t48;
t68 = -t47 * t25 + t49 * t29;
t66 = t48 * mrSges(5,1) - t46 * mrSges(5,2);
t50 = sin(qJ(5));
t52 = cos(qJ(5));
t65 = -t46 * t52 - t48 * t50;
t64 = t46 * t50 - t48 * t52;
t62 = -qJDD(1) * mrSges(5,3) - t54 * t66;
t10 = (pkin(6) * qJDD(1) + t48 * t77 - t16) * t46 + t74;
t11 = -pkin(6) * t72 - t41 * t77 + t70;
t31 = t64 * qJD(1);
t32 = t65 * qJD(1);
t18 = -t31 * mrSges(6,1) + t32 * mrSges(6,2);
t21 = t31 * qJD(5) + t65 * qJDD(1);
t26 = -qJD(5) * mrSges(6,2) + t31 * mrSges(6,3);
t8 = m(6) * (t52 * t10 - t50 * t11) - t21 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t32 * t18 + qJD(5) * t26;
t20 = -t32 * qJD(5) + t64 * qJDD(1);
t27 = qJD(5) * mrSges(6,1) - t32 * mrSges(6,3);
t9 = m(6) * (t50 * t10 + t52 * t11) + t20 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t31 * t18 - qJD(5) * t27;
t5 = m(5) * t74 + t50 * t9 + t52 * t8 + (-m(5) * t16 - t62) * t46;
t6 = m(5) * t70 + t62 * t48 - t50 * t8 + t52 * t9;
t4 = m(4) * t75 - t54 * mrSges(4,1) + qJDD(1) * mrSges(4,2) - t46 * t5 + t48 * t6;
t61 = qJDD(1) * pkin(3) + qJDD(4) - t68;
t56 = t20 * mrSges(6,1) + t31 * t26 - m(6) * (pkin(4) * t72 + (-t73 * pkin(6) - qJ(4)) * t54 + t61) - t32 * t27 - t21 * mrSges(6,2);
t55 = m(5) * (-t54 * qJ(4) + t61) - t56;
t7 = m(4) * t68 + (-mrSges(4,2) + t80) * t54 + (-mrSges(4,1) - t66) * qJDD(1) - t55;
t63 = t49 * t4 - t47 * t7 + m(3) * (-t54 * pkin(1) + t60) + qJDD(1) * mrSges(3,3);
t59 = -m(3) * (-qJDD(1) * pkin(1) + t57) - t47 * t4 - t49 * t7;
t58 = m(4) * t44 + t46 * t6 + t48 * t5;
t2 = m(2) * t69 + (-mrSges(2,2) + mrSges(3,3)) * t54 + t76 * qJDD(1) + t59;
t1 = m(2) * t67 - qJDD(1) * mrSges(2,2) - t76 * t54 + t63;
t3 = [-m(1) * g(1) + t53 * t1 - t51 * t2, t1, -t54 * mrSges(3,1) + t63, t4, t6, t9; -m(1) * g(2) + t51 * t1 + t53 * t2, t2, -m(3) * g(3) - t58, t7, t5, t8; (-m(1) + t79) * g(3) - t58, t79 * g(3) - t58, -qJDD(1) * mrSges(3,1) - t54 * mrSges(3,3) - t59, t58, t66 * qJDD(1) - t54 * t80 + t55, -t56;];
f_new = t3;
