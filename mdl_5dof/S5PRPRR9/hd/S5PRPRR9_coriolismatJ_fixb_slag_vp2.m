% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRR9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR9_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR9_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:39
% EndTime: 2019-12-31 17:39:39
% DurationCPUTime: 0.32s
% Computational Cost: add. (769->90), mult. (1443->125), div. (0->0), fcn. (920->4), ass. (0->58)
t59 = cos(qJ(5));
t96 = -t59 / 0.2e1;
t60 = cos(qJ(4));
t57 = sin(qJ(5));
t82 = t57 ^ 2 + t59 ^ 2;
t73 = t82 * t60;
t75 = Ifges(6,1) * t59 / 0.2e1 - Ifges(6,4) * t57 + Ifges(6,2) * t96;
t95 = t75 * t57;
t44 = -t59 * mrSges(6,1) + t57 * mrSges(6,2);
t58 = sin(qJ(4));
t37 = t58 * t44;
t72 = t60 * mrSges(5,2) - mrSges(6,3) * t73;
t94 = -t58 * mrSges(5,1) + t37 - t72;
t61 = -pkin(2) - pkin(3);
t42 = -t58 * qJ(3) + t60 * t61;
t93 = -t42 / 0.2e1;
t91 = m(6) * (-0.1e1 + t82) * t60 * t58;
t84 = t59 * mrSges(6,2);
t85 = t57 * mrSges(6,1);
t45 = t84 + t85;
t90 = pkin(4) * t45;
t87 = Ifges(6,4) * t59;
t40 = pkin(4) - t42;
t86 = t40 * t45;
t83 = t60 * t45;
t43 = t60 * qJ(3) + t58 * t61;
t41 = -pkin(7) + t43;
t74 = t82 * t42;
t63 = mrSges(6,3) * t74 - t42 * mrSges(5,2) + (-mrSges(5,1) + t44) * t43;
t2 = -m(6) * (t40 * t43 + t41 * t74) + t63;
t81 = t2 * qJD(2);
t6 = -mrSges(4,3) - m(4) * qJ(3) - m(6) * (t40 * t58 + t41 * t73) - m(5) * (-t42 * t58 + t43 * t60) + t94;
t80 = t6 * qJD(2);
t47 = -Ifges(6,2) * t57 + t87;
t48 = Ifges(6,1) * t57 + t87;
t76 = t48 / 0.2e1 + t47 / 0.2e1;
t64 = t76 * t59 + t95;
t9 = t64 - t86;
t79 = t9 * qJD(2);
t70 = -t84 / 0.2e1 - t85 / 0.2e1;
t16 = (-t45 / 0.2e1 + t70) * t60;
t78 = t16 * qJD(2);
t77 = qJD(3) * t91;
t71 = -Ifges(6,5) * t59 + Ifges(6,6) * t57;
t62 = -t37 / 0.2e1 + m(6) * ((t82 * t41 - t43) * t60 + (t40 + t74) * t58) / 0.2e1;
t65 = m(6) * (-pkin(4) * t58 + pkin(7) * t73);
t1 = (mrSges(5,1) - t44 / 0.2e1) * t58 - t65 / 0.2e1 + t62 + t72;
t69 = t1 * qJD(2) + t77;
t67 = t70 * t60;
t10 = t64 - t90;
t15 = (t45 / 0.2e1 + t70) * t60;
t4 = (-pkin(4) / 0.2e1 - t40 / 0.2e1) * t45 + (mrSges(6,2) * t93 + t76) * t59 + (mrSges(6,1) * t93 + t75) * t57;
t66 = t4 * qJD(2) + t15 * qJD(3) - t10 * qJD(4);
t18 = t83 / 0.2e1 + t67;
t17 = -t83 / 0.2e1 + t67;
t5 = t90 / 0.2e1 + t86 / 0.2e1 + t70 * t42 + (t48 + t47) * t96 - t95;
t3 = t65 / 0.2e1 + t37 / 0.2e1 + t62;
t7 = [0, 0, 0, 0, t45 * qJD(5); 0, -t6 * qJD(3) - t2 * qJD(4) + t9 * qJD(5), t3 * qJD(4) + t18 * qJD(5) + t77 - t80, -t81 + t3 * qJD(3) + (m(6) * (-pkin(4) * t43 + pkin(7) * t74) + t63) * qJD(4) + t5 * qJD(5), t79 + t18 * qJD(3) + t5 * qJD(4) + (t44 * t41 + t71) * qJD(5); 0, t1 * qJD(4) - t16 * qJD(5) + t80, qJD(4) * t91, (t65 + t94) * qJD(4) + t17 * qJD(5) + t69, t17 * qJD(4) + qJD(5) * t37 - t78; 0, -t1 * qJD(3) - t4 * qJD(5) + t81, -t15 * qJD(5) - t69, t10 * qJD(5), (t44 * pkin(7) - t71) * qJD(5) - t66; 0, t16 * qJD(3) + t4 * qJD(4) - t79, t15 * qJD(4) + t78, t66, 0;];
Cq = t7;
