% Return the minimum parameter vector for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% MPV [19x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPPPR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t83 = sin(pkin(7));
t86 = cos(pkin(7));
t93 = t86 * mrSges(3,1) - t83 * mrSges(3,2);
t90 = 2 * pkin(6) * mrSges(6,3) + Ifges(6,2);
t89 = -pkin(6) * m(6) - mrSges(6,3);
t88 = (pkin(4) ^ 2);
t87 = pkin(6) ^ 2;
t85 = cos(pkin(8));
t84 = cos(pkin(9));
t82 = sin(pkin(8));
t81 = sin(pkin(9));
t1 = [Ifges(2,3) + Ifges(3,3) + t85 ^ 2 * (t88 * m(6) + Ifges(4,2) + Ifges(5,3)) + (0.2e1 * t85 * (Ifges(4,4) - t84 * (t89 * pkin(4) + Ifges(5,5)) + t81 * Ifges(5,6)) + (Ifges(4,1) + t84 ^ 2 * (m(6) * t87 + Ifges(5,1) + t90) + (-0.2e1 * t84 * Ifges(5,4) + (Ifges(5,2) + (t87 + t88) * m(6) + t90) * t81) * t81) * t82) * t82 + 0.2e1 * t93 * pkin(1); mrSges(2,1) + t93; t83 * mrSges(3,1) + t86 * mrSges(3,2) + mrSges(2,2); m(3); mrSges(4,1); mrSges(4,2); mrSges(4,3); m(4); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t89; m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
