% Return the minimum parameter vector for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% MPV [22x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5PRPRR6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR6_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t110 = (m(5) + m(6));
t114 = (pkin(4) ^ 2);
t118 = (t114 * m(6) + Ifges(5,2));
t117 = 2 * pkin(8) * mrSges(6,3) + Ifges(6,2);
t116 = 2 * pkin(7) * mrSges(5,3) + t118;
t115 = pkin(8) * m(6) + mrSges(6,3);
t113 = pkin(7) ^ 2;
t112 = pkin(8) ^ 2;
t109 = cos(pkin(10));
t108 = sin(pkin(10));
t1 = [m(2) + m(3); Ifges(3,3) + t109 ^ 2 * (Ifges(4,2) + (pkin(3) ^ 2 + t113) * t110 + t116) + (0.2e1 * t109 * Ifges(4,4) + (t113 * t110 + Ifges(4,1) + t116) * t108) * t108; mrSges(3,1); mrSges(3,2); pkin(3) * t110 + mrSges(4,1); mrSges(4,2); pkin(7) * t110 + mrSges(4,3) + mrSges(5,3); m(4) + t110; m(6) * t112 + Ifges(5,1) + t117 - t118; t115 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t112 + t114) * m(6) + t117; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t115; Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
