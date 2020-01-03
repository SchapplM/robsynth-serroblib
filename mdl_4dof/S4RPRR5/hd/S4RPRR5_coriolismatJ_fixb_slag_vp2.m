% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:33
% EndTime: 2019-12-31 16:51:34
% DurationCPUTime: 0.32s
% Computational Cost: add. (768->91), mult. (1440->124), div. (0->0), fcn. (918->4), ass. (0->62)
t59 = cos(qJ(4));
t101 = -t59 / 0.2e1;
t57 = sin(qJ(4));
t76 = Ifges(5,1) * t59 / 0.2e1 - Ifges(5,4) * t57 + Ifges(5,2) * t101;
t100 = t76 * t57;
t84 = t57 ^ 2 + t59 ^ 2;
t99 = mrSges(5,3) * t84;
t58 = sin(qJ(3));
t82 = t58 * qJ(2);
t60 = cos(qJ(3));
t61 = -pkin(1) - pkin(2);
t86 = t60 * t61;
t42 = -t82 + t86;
t98 = -t42 / 0.2e1;
t96 = m(5) * (-0.1e1 + t84) * t60 * t58;
t74 = t84 * t60;
t95 = m(5) * (-pkin(3) * t58 + pkin(6) * t74);
t89 = t59 * mrSges(5,2);
t90 = t57 * mrSges(5,1);
t45 = t89 + t90;
t94 = pkin(3) * t45;
t92 = Ifges(5,4) * t59;
t44 = -t59 * mrSges(5,1) + t57 * mrSges(5,2);
t38 = t58 * t44;
t87 = t60 * t45;
t85 = -t58 * mrSges(4,1) - t60 * mrSges(4,2);
t43 = t60 * qJ(2) + t58 * t61;
t41 = -pkin(6) + t43;
t71 = pkin(3) - t42;
t72 = -t42 * mrSges(4,2) + (-mrSges(4,1) + t44) * t43;
t75 = t84 * t42;
t2 = -m(5) * (t41 * t75 + t71 * t43) + t72 + t42 * t99;
t83 = t2 * qJD(1);
t64 = t60 * t99 + t85;
t6 = t38 - m(5) * (t41 * t74 + t71 * t58) - m(4) * (-t42 * t58 + t43 * t60) - mrSges(3,3) - m(3) * qJ(2) + t64;
t81 = t6 * qJD(1);
t47 = -Ifges(5,2) * t57 + t92;
t48 = Ifges(5,1) * t57 + t92;
t77 = t48 / 0.2e1 + t47 / 0.2e1;
t63 = t77 * t59 + t100;
t67 = t71 * t45;
t9 = -t67 + t63;
t80 = t9 * qJD(1);
t69 = -t89 / 0.2e1 - t90 / 0.2e1;
t16 = (-t45 / 0.2e1 + t69) * t60;
t79 = t16 * qJD(1);
t78 = qJD(2) * t96;
t73 = -Ifges(5,5) * t59 + Ifges(5,6) * t57;
t70 = t95 / 0.2e1 + t38 / 0.2e1;
t62 = -t38 / 0.2e1 + m(5) * ((t84 * t41 - t43) * t60 + (t75 + t71) * t58) / 0.2e1;
t1 = -t62 + t64 + t70;
t68 = -t1 * qJD(1) + t78;
t66 = t69 * t60;
t10 = t63 - t94;
t15 = (t45 / 0.2e1 + t69) * t60;
t4 = (-pkin(3) - t82 / 0.2e1 + t86 / 0.2e1) * t45 + (mrSges(5,2) * t98 + t77) * t59 + (mrSges(5,1) * t98 + t76) * t57;
t65 = t4 * qJD(1) + t15 * qJD(2) - t10 * qJD(3);
t18 = t87 / 0.2e1 + t66;
t17 = -t87 / 0.2e1 + t66;
t5 = t94 / 0.2e1 + t67 / 0.2e1 + t69 * t42 + (t48 + t47) * t101 - t100;
t3 = t62 + t70;
t7 = [-t6 * qJD(2) - t2 * qJD(3) + t9 * qJD(4), t3 * qJD(3) + t18 * qJD(4) + t78 - t81, -t83 + t3 * qJD(2) + (m(5) * (-pkin(3) * t43 + pkin(6) * t75) + mrSges(5,3) * t75 + t72) * qJD(3) + t5 * qJD(4), t80 + t18 * qJD(2) + t5 * qJD(3) + (t44 * t41 + t73) * qJD(4); -t1 * qJD(3) - t16 * qJD(4) + t81, qJD(3) * t96, (mrSges(5,3) * t74 + t38 + t85 + t95) * qJD(3) + t17 * qJD(4) + t68, t17 * qJD(3) + qJD(4) * t38 - t79; t1 * qJD(2) - t4 * qJD(4) + t83, -t15 * qJD(4) - t68, t10 * qJD(4), (t44 * pkin(6) - t73) * qJD(4) - t65; t16 * qJD(2) + t4 * qJD(3) - t80, t15 * qJD(3) + t79, t65, 0;];
Cq = t7;
