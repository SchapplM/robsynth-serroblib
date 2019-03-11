% Return the minimum parameter vector for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [35x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRRR12_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_convert_par2_MPV_fixb: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR12_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR12_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t247 = (m(6) + m(7));
t236 = m(5) + t247;
t229 = m(4) + t236;
t278 = pkin(10) * t229;
t241 = sin(pkin(7));
t244 = cos(pkin(8));
t262 = pkin(11) * t236 + mrSges(5,3);
t274 = t262 * t244;
t223 = mrSges(4,3) + t274;
t260 = t223 + t278;
t277 = t260 * t241;
t245 = cos(pkin(7));
t276 = t260 * t245;
t240 = sin(pkin(8));
t275 = t262 * t240;
t234 = t244 ^ 2;
t254 = (pkin(4) ^ 2);
t228 = t254 * t247 + Ifges(5,2);
t251 = pkin(11) ^ 2;
t255 = pkin(3) ^ 2;
t265 = Ifges(4,2) + (t251 * t234 + t255) * t236 + t234 * t228;
t272 = (pkin(11) * mrSges(5,3));
t267 = 2 * t272;
t257 = t234 * t267 + (0.2e1 * t223 + t278) * pkin(10) + t265;
t273 = t251 * t236 + t228;
t253 = (pkin(5) ^ 2);
t271 = (t253 * m(7) + Ifges(6,2));
t270 = 2 * pkin(13) * mrSges(7,3) + Ifges(7,2);
t268 = pkin(2) ^ 2 * t229;
t264 = 2 * pkin(12) * mrSges(6,3) + t271;
t263 = pkin(13) * m(7) + mrSges(7,3);
t261 = pkin(12) * t247 + mrSges(6,3);
t258 = t267 + t273;
t250 = pkin(12) ^ 2;
t249 = pkin(13) ^ 2;
t246 = cos(pkin(6));
t243 = cos(pkin(14));
t242 = sin(pkin(6));
t239 = sin(pkin(14));
t1 = [Ifges(2,3) + t246 ^ 2 * (t257 * t241 ^ 2 + Ifges(3,3) + t268) + (0.2e1 * (t239 * (-pkin(2) * t276 + Ifges(3,5)) + t243 * (t257 * t245 * t241 + Ifges(3,6))) * t246 + (t243 ^ 2 * (t257 * t245 ^ 2 + Ifges(3,2) + t268) + (0.2e1 * t243 * (pkin(2) * t277 + Ifges(3,4)) + (Ifges(3,1) + t257) * t239) * t239) * t242) * t242; mrSges(2,1); mrSges(2,2); pkin(2) * t229 + mrSges(3,1); mrSges(3,2) - t277; mrSges(3,3) + t276; m(3) + t229; Ifges(4,1) + (-0.2e1 * t234 + 0.2e1) * t272 - t265 + t273; pkin(3) * t275 + Ifges(4,4); -pkin(3) * t274 + Ifges(4,5); t258 * t244 * t240 + Ifges(4,6); t258 * t240 ^ 2 + t255 * t236 + Ifges(4,3); pkin(3) * t236 + mrSges(4,1); mrSges(4,2) - t275; (t250 * t247) + Ifges(5,1) - t228 + t264; t261 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t250 + t254) * t247 + t264; pkin(4) * t247 + mrSges(5,1); mrSges(5,2) - t261; m(7) * t249 + Ifges(6,1) + t270 - t271; t263 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t249 + t253) * m(7) + t270; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t263; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
