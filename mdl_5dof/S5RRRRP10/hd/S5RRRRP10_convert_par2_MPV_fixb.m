% Return the minimum parameter vector for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% MPV [26x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRRP10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP10_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t123 = (m(4) + m(5));
t119 = m(3) + t123;
t135 = (t119 * pkin(7));
t134 = (-Ifges(5,2) - Ifges(6,2));
t127 = (pkin(3) ^ 2);
t133 = (t127 * m(5) + Ifges(4,2));
t132 = 2 * pkin(9) * mrSges(5,3) - t134;
t131 = 2 * pkin(8) * mrSges(4,3) + t133;
t130 = pkin(9) * m(5) + mrSges(5,3);
t129 = pkin(8) * t123 + mrSges(4,3);
t128 = (pkin(2) ^ 2);
t126 = pkin(8) ^ 2;
t125 = pkin(9) ^ 2;
t122 = sin(pkin(5));
t1 = [pkin(1) ^ 2 * t119 + Ifges(2,3) + (t128 * t123 + Ifges(3,2) + (2 * mrSges(3,3) + t135) * pkin(7)) * t122 ^ 2; pkin(1) * t119 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t135) * t122; Ifges(3,1) - Ifges(3,2) + (t126 - t128) * t123 + t131; t129 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + (t126 + t128) * t123 + t131; pkin(2) * t123 + mrSges(3,1); mrSges(3,2) - t129; t125 * m(5) + Ifges(4,1) + t132 - t133; t130 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t125 + t127) * m(5) + t132; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t130; Ifges(5,1) + Ifges(6,1) + t134; Ifges(5,4) + Ifges(6,4); Ifges(5,5) + Ifges(6,5); Ifges(5,6) + Ifges(6,6); 2 * pkin(4) * mrSges(6,1) + Ifges(5,3) + Ifges(6,3); mrSges(5,1) + mrSges(6,1); mrSges(5,2) + mrSges(6,2); mrSges(6,3); m(6);];
MPV = t1;
