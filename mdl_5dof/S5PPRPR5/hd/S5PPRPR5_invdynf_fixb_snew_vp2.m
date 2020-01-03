% Calculate vector of cutting forces with Newton-Euler
% S5PPRPR5
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPRPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:25
% EndTime: 2019-12-31 17:33:26
% DurationCPUTime: 0.20s
% Computational Cost: add. (1134->65), mult. (1817->85), div. (0->0), fcn. (908->6), ass. (0->43)
t59 = -pkin(3) - pkin(6);
t58 = mrSges(4,1) - mrSges(5,2);
t57 = -mrSges(4,2) + mrSges(5,3);
t33 = sin(pkin(7));
t34 = cos(pkin(7));
t49 = t33 * g(1) - t34 * g(2);
t24 = qJDD(2) - t49;
t26 = -t34 * g(1) - t33 * g(2);
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t56 = t36 * t24 + t38 * t26;
t35 = sin(qJ(5));
t55 = qJD(3) * t35;
t37 = cos(qJ(5));
t54 = qJD(3) * t37;
t53 = qJD(3) * qJD(5);
t32 = g(3) - qJDD(1);
t39 = qJD(3) ^ 2;
t51 = t38 * t24 - t36 * t26;
t42 = -t39 * qJ(4) + qJDD(4) - t51;
t12 = t59 * qJDD(3) + t42;
t21 = (mrSges(6,1) * t35 + mrSges(6,2) * t37) * qJD(3);
t23 = t37 * qJDD(3) - t35 * t53;
t27 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t55;
t8 = m(6) * (t37 * t12 - t35 * t32) - t23 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t21 * t54 + qJD(5) * t27;
t22 = -t35 * qJDD(3) - t37 * t53;
t28 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t54;
t9 = m(6) * (t35 * t12 + t37 * t32) + t22 * mrSges(6,3) - qJDD(5) * mrSges(6,2) - t21 * t55 - qJD(5) * t28;
t52 = m(5) * t32 - t35 * t8 + t37 * t9;
t43 = -m(5) * (-qJDD(3) * pkin(3) + t42) - t35 * t9 - t37 * t8;
t3 = m(4) * t51 + t58 * qJDD(3) + t57 * t39 + t43;
t41 = qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t56;
t46 = -t22 * mrSges(6,1) + m(6) * (t59 * t39 + t41) + t27 * t55 + t28 * t54 + t23 * mrSges(6,2);
t40 = -m(5) * (t39 * pkin(3) - t41) + t46;
t5 = m(4) * t56 + t57 * qJDD(3) - t58 * t39 + t40;
t50 = m(3) * t26 - t36 * t3 + t38 * t5;
t48 = m(4) * t32 + t52;
t47 = -m(3) * t32 - t48;
t45 = -m(2) * t32 + t47;
t44 = m(3) * t24 + t38 * t3 + t36 * t5;
t2 = m(2) * t26 + t50;
t1 = m(2) * t49 - t44;
t4 = [-m(1) * g(1) - t33 * t1 + t34 * t2, t2, t50, t5, t52, t9; -m(1) * g(2) + t34 * t1 + t33 * t2, t1, t47, t3, -t39 * mrSges(5,2) - qJDD(3) * mrSges(5,3) - t40, t8; -m(1) * g(3) + t45, t45, t44, t48, qJDD(3) * mrSges(5,2) - t39 * mrSges(5,3) - t43, t46;];
f_new = t4;
