% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:32
% EndTime: 2019-07-18 18:16:32
% DurationCPUTime: 0.27s
% Computational Cost: add. (857->71), mult. (1345->89), div. (0->0), fcn. (863->4), ass. (0->53)
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t67 = t46 * mrSges(5,1) + t48 * mrSges(5,2);
t83 = qJD(3) * t67;
t82 = t67 * qJD(4);
t47 = sin(qJ(2));
t42 = t47 * pkin(1) + qJ(3);
t81 = m(4) * t42;
t49 = cos(qJ(2));
t80 = -t49 * mrSges(3,2) + (-mrSges(3,1) - mrSges(4,1)) * t47;
t32 = (t46 * t47 + t48 * t49) * pkin(1);
t30 = t32 * mrSges(5,2);
t77 = -t30 / 0.2e1;
t75 = t49 * pkin(1);
t74 = m(4) * qJ(3);
t60 = -pkin(2) - t75;
t40 = -pkin(3) + t60;
t27 = t46 * t40 + t48 * t42;
t73 = t27 * mrSges(5,1);
t22 = t27 * t48;
t31 = (-t46 * t49 + t47 * t48) * pkin(1);
t72 = t31 * mrSges(5,1);
t50 = -pkin(2) - pkin(3);
t37 = -t46 * qJ(3) + t48 * t50;
t71 = t37 * mrSges(5,2);
t38 = t48 * qJ(3) + t46 * t50;
t33 = t38 * t48;
t26 = t48 * t40 - t46 * t42;
t57 = mrSges(4,3) * t75 + t30 - t72;
t3 = -m(5) * (t26 * t31 + t27 * t32) - t57 + (-m(4) * (t42 * t49 + t60 * t47) - t80) * pkin(1);
t66 = t3 * qJD(1);
t25 = t26 * mrSges(5,2);
t7 = -t25 - t73;
t65 = t7 * qJD(1);
t64 = t7 * qJD(4);
t36 = t38 * mrSges(5,1);
t10 = t36 + t71;
t63 = t10 * qJD(4);
t61 = mrSges(4,3) + t67;
t12 = m(5) * (-t26 * t46 + t22) + t81 + t61;
t62 = t12 * qJD(1);
t53 = -t25 / 0.2e1 - t36 / 0.2e1 - t71 / 0.2e1;
t1 = t77 + (t31 / 0.2e1 - t27 / 0.2e1) * mrSges(5,1) + t53;
t56 = -t1 * qJD(1) + t10 * qJD(2);
t16 = t74 + m(5) * (-t37 * t46 + t33) + t61;
t51 = -m(5) * (t22 + t33 + (-t26 - t37) * t46) / 0.2e1 - t61;
t52 = m(5) * (t48 * t31 + t46 * t32) / 0.2e1;
t4 = t51 + t52 - t74;
t55 = -t4 * qJD(1) + t16 * qJD(2);
t54 = (-qJD(1) - qJD(2)) * t67;
t5 = t52 - t51 + t81;
t2 = t73 / 0.2e1 + t72 / 0.2e1 + t77 - t53;
t6 = [-t3 * qJD(2) + t12 * qJD(3) - t64, t5 * qJD(3) + t2 * qJD(4) - t66 + (m(5) * (t37 * t31 + t38 * t32) + t57 + (m(4) * (-pkin(2) * t47 + qJ(3) * t49) + t80) * pkin(1)) * qJD(2), t5 * qJD(2) + t62, t2 * qJD(2) + t64 - t65; -t4 * qJD(3) - t1 * qJD(4) + t66, t16 * qJD(3) + t63, t55, t56 - t63; t4 * qJD(2) - t62 + t82, -t55 + t82, 0, -t54 - t82; t1 * qJD(2) + t65 - t83, -t56 - t83, t54, 0;];
Cq  = t6;
