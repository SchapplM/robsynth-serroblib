% Return the minimum parameter vector for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRPR10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR10_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t139 = m(3) + m(4);
t152 = (t139 * pkin(7));
t151 = 2 * pkin(9) * mrSges(6,3) + Ifges(6,2);
t136 = sin(pkin(10));
t138 = cos(pkin(10));
t150 = t136 * t138;
t140 = pkin(9) ^ 2;
t125 = m(6) * t140 + Ifges(5,1) + t151;
t142 = (pkin(4) ^ 2);
t129 = t142 * m(6) + Ifges(5,2);
t131 = t136 ^ 2;
t133 = t138 ^ 2;
t149 = t131 * t125 + t133 * t129 + Ifges(4,2);
t148 = pkin(8) * m(4) + mrSges(4,3);
t147 = pkin(9) * m(6) + mrSges(6,3);
t127 = t147 * pkin(4) + Ifges(5,4);
t146 = t127 * t150;
t145 = (2 * pkin(8) * mrSges(4,3)) + 0.2e1 * t146 + t149;
t128 = mrSges(5,2) - t147;
t130 = m(6) * pkin(4) + mrSges(5,1);
t144 = -t136 * t128 + t138 * t130;
t143 = (pkin(2) ^ 2);
t141 = pkin(8) ^ 2;
t137 = sin(pkin(5));
t1 = [pkin(1) ^ 2 * t139 + Ifges(2,3) + (t143 * m(4) + Ifges(3,2) + (2 * mrSges(3,3) + t152) * pkin(7)) * t137 ^ 2; pkin(1) * t139 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t152) * t137; Ifges(3,1) - Ifges(3,2) + ((t141 - t143) * m(4)) + t145; t148 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + ((t141 + t143) * m(4)) + t145; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t148; t133 * t125 + t131 * t129 + Ifges(4,1) - 0.4e1 * t146 - t149; Ifges(4,4) + (t133 - t131) * t127 + (t125 - t129) * t150; t138 * Ifges(5,5) - t136 * Ifges(5,6) + Ifges(4,5); t136 * Ifges(5,5) + t138 * Ifges(5,6) + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + ((t140 + t142) * m(6)) + 0.2e1 * t144 * pkin(3) + t151; mrSges(4,1) + t144; t138 * t128 + t136 * t130 + mrSges(4,2); mrSges(5,3); m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
