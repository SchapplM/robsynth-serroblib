% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPPR2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(7)) + r_i_i_C(2) * sin(pkin(7)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:31
	% EndTime: 2022-01-23 09:00:31
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (21->14), mult. (54->21), div. (0->0), fcn. (44->6), ass. (0->11)
	t118 = sin(qJ(1));
	t124 = qJD(1) * t118;
	t119 = cos(qJ(1));
	t123 = qJD(1) * t119;
	t115 = sin(pkin(7));
	t122 = t115 * qJD(3);
	t114 = sin(pkin(8));
	t116 = cos(pkin(8));
	t121 = t114 * r_i_i_C(1) + t116 * r_i_i_C(2) + qJ(2);
	t120 = -pkin(1) + (-r_i_i_C(1) * t116 + r_i_i_C(2) * t114 - pkin(2)) * cos(pkin(7)) + (-r_i_i_C(3) - qJ(3)) * t115;
	t1 = [-t118 * t122 + t119 * qJD(2) + (-t121 * t118 + t120 * t119) * qJD(1), t123, -t115 * t124, 0, 0; t119 * t122 + t118 * qJD(2) + (t120 * t118 + t121 * t119) * qJD(1), t124, t115 * t123, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:32
	% EndTime: 2022-01-23 09:00:32
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (41->26), mult. (100->45), div. (0->0), fcn. (90->8), ass. (0->20)
	t172 = cos(pkin(7));
	t173 = sin(qJ(1));
	t180 = t172 * t173;
	t174 = cos(qJ(1));
	t179 = t172 * t174;
	t178 = qJD(1) * t173;
	t177 = qJD(1) * t174;
	t167 = sin(pkin(9));
	t170 = cos(pkin(9));
	t176 = qJD(1) * (t170 * r_i_i_C(1) - t167 * r_i_i_C(2));
	t168 = sin(pkin(8));
	t169 = sin(pkin(7));
	t171 = cos(pkin(8));
	t175 = -(t171 * pkin(3) + t168 * qJ(4) + pkin(2)) * t172 - pkin(1) + (-r_i_i_C(1) * t167 - r_i_i_C(2) * t170 - qJ(3)) * t169;
	t166 = qJD(4) * t171 - qJD(2);
	t165 = -t168 * pkin(3) + qJ(4) * t171 - qJ(2);
	t164 = t168 * qJD(4) * t172 + t169 * qJD(3);
	t162 = (t168 * t179 - t171 * t173) * qJD(1);
	t160 = (-t168 * t180 - t171 * t174) * qJD(1);
	t1 = [-t162 * r_i_i_C(3) - t164 * t173 - t166 * t174 + (-t168 * t173 - t171 * t179) * t176 + (t165 * t173 + t175 * t174) * qJD(1), t177, -t169 * t178, t160, 0; t160 * r_i_i_C(3) + t164 * t174 - t166 * t173 + (t168 * t174 - t171 * t180) * t176 + (-t165 * t174 + t175 * t173) * qJD(1), t178, t169 * t177, t162, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:00:32
	% EndTime: 2022-01-23 09:00:32
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (126->64), mult. (332->115), div. (0->0), fcn. (350->10), ass. (0->47)
	t276 = sin(pkin(9));
	t281 = cos(pkin(7));
	t319 = t281 * t276;
	t279 = cos(pkin(9));
	t273 = t279 * pkin(4) + t276 * pkin(6) + pkin(3);
	t277 = sin(pkin(8));
	t278 = sin(pkin(7));
	t280 = cos(pkin(8));
	t318 = -qJD(4) * t280 + qJD(2) - ((t277 * qJ(4) + t273 * t280 + pkin(2)) * t281 + pkin(1) + (t276 * pkin(4) - t279 * pkin(6) + qJ(3)) * t278) * qJD(1);
	t302 = t279 * t280;
	t269 = t278 * t276 + t281 * t302;
	t282 = sin(qJ(5));
	t284 = cos(qJ(5));
	t304 = t277 * t284;
	t313 = -t269 * t282 + t281 * t304;
	t259 = t313 * qJD(5);
	t305 = t277 * t282;
	t270 = t279 * t305 + t284 * t280;
	t264 = t270 * qJD(5);
	t283 = sin(qJ(1));
	t285 = cos(qJ(1));
	t317 = -t259 * t283 - t285 * t264;
	t271 = t279 * t304 - t282 * t280;
	t290 = t269 * t284 + t281 * t305;
	t316 = -t271 * t283 - t285 * t290;
	t300 = t279 * t285;
	t314 = t284 * (-t269 * t283 + t277 * t300);
	t312 = -t285 * t270 - t283 * t313;
	t260 = t290 * qJD(5);
	t265 = t271 * qJD(5);
	t311 = -t260 * t285 - t283 * t265;
	t310 = r_i_i_C(3) * qJD(1);
	t303 = t277 * t285;
	t301 = t279 * t283;
	t299 = t281 * t283;
	t296 = t285 * t280;
	t295 = qJD(1) * t283;
	t294 = qJD(1) * t285;
	t289 = -t280 * t283 + t281 * t303;
	t288 = t277 * t299 + t296;
	t287 = qJ(4) * t280 - t273 * t277 - qJ(2);
	t286 = (t269 * t285 + t277 * t301) * t282;
	t272 = t277 * qJD(4) * t281 + qJD(3) * t278;
	t268 = t278 * t302 - t319;
	t267 = t289 * qJD(1);
	t266 = t288 * qJD(1);
	t1 = [(t316 * qJD(1) + t317) * r_i_i_C(1) + (t278 * t300 - t296 * t319) * t310 + t287 * t295 + ((t270 * qJD(1) + t260) * r_i_i_C(2) - t277 * t276 * t310 - t272) * t283 + ((-t313 * qJD(1) - t265) * r_i_i_C(2) + t318) * t285, t294, -t278 * t295, -t266, t311 * r_i_i_C(1) + (-t259 * t285 + t283 * t264) * r_i_i_C(2) + (t312 * r_i_i_C(1) + (-t271 * t285 + t283 * t290) * r_i_i_C(2)) * qJD(1); (-t266 * t282 + (t284 * t289 - t286) * qJD(5) + qJD(1) * t314) * r_i_i_C(1) + (t312 * qJD(1) + t311) * r_i_i_C(2) + ((-t280 * t299 + t303) * t276 + t278 * t301) * t310 + t272 * t285 - t287 * t294 + t318 * t283, t295, t278 * t294, t267, (t267 * t284 + (-t282 * t288 + t314) * qJD(5)) * r_i_i_C(1) + t317 * r_i_i_C(2) + (-r_i_i_C(1) * t286 + t316 * r_i_i_C(2)) * qJD(1); 0, 0, 0, 0, ((-t268 * t284 - t278 * t305) * r_i_i_C(1) + (t268 * t282 - t278 * t304) * r_i_i_C(2)) * qJD(5);];
	JaD_transl = t1;
end