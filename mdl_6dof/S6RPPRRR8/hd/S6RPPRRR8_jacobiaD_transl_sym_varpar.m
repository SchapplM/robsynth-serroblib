% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRR8
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:13
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->6), mult. (20->10), div. (0->0), fcn. (12->2), ass. (0->5)
	t10 = -pkin(1) + r_i_i_C(2);
	t9 = r_i_i_C(3) + qJ(2);
	t8 = cos(qJ(1));
	t7 = sin(qJ(1));
	t1 = [t8 * qJD(2) + (t10 * t8 - t9 * t7) * qJD(1), qJD(1) * t8, 0, 0, 0, 0; t7 * qJD(2) + (t10 * t7 + t9 * t8) * qJD(1), qJD(1) * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (15->10), mult. (36->14), div. (0->0), fcn. (24->4), ass. (0->7)
	t14 = sin(qJ(1));
	t18 = qJD(1) * t14;
	t17 = -pkin(1) - r_i_i_C(3) - qJ(3);
	t16 = r_i_i_C(1) * sin(pkin(10)) + r_i_i_C(2) * cos(pkin(10)) + qJ(2);
	t15 = cos(qJ(1));
	t11 = qJD(1) * t15;
	t1 = [t15 * qJD(2) - t14 * qJD(3) + (-t16 * t14 + t17 * t15) * qJD(1), t11, -t18, 0, 0, 0; t14 * qJD(2) + t15 * qJD(3) + (t17 * t14 + t16 * t15) * qJD(1), t18, t11, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (48->21), mult. (82->32), div. (0->0), fcn. (54->5), ass. (0->14)
	t24 = pkin(10) + qJ(4);
	t21 = sin(t24);
	t22 = cos(t24);
	t36 = (r_i_i_C(1) * t22 - r_i_i_C(2) * t21) * qJD(4);
	t27 = sin(qJ(1));
	t35 = qJD(1) * t27;
	t28 = cos(qJ(1));
	t23 = qJD(1) * t28;
	t34 = qJD(4) * t27;
	t33 = qJD(4) * t28;
	t32 = -pkin(1) - r_i_i_C(3) - pkin(7) - qJ(3);
	t30 = pkin(3) * sin(pkin(10)) + r_i_i_C(1) * t21 + r_i_i_C(2) * t22 + qJ(2);
	t29 = qJD(2) + t36;
	t1 = [-t27 * qJD(3) + t29 * t28 + (-t30 * t27 + t32 * t28) * qJD(1), t23, -t35, (-t21 * t23 - t22 * t34) * r_i_i_C(2) + (-t21 * t34 + t22 * t23) * r_i_i_C(1), 0, 0; t28 * qJD(3) + t29 * t27 + (t32 * t27 + t30 * t28) * qJD(1), t35, t23, (-t21 * t35 + t22 * t33) * r_i_i_C(2) + (t21 * t33 + t22 * t35) * r_i_i_C(1), 0, 0; 0, 0, 0, -t36, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:29
	% EndTime: 2019-10-10 00:13:30
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (165->45), mult. (288->78), div. (0->0), fcn. (225->7), ass. (0->33)
	t211 = sin(qJ(5));
	t213 = cos(qJ(5));
	t217 = (r_i_i_C(1) * t211 + r_i_i_C(2) * t213) * qJD(5);
	t232 = pkin(8) + r_i_i_C(3);
	t238 = t232 * qJD(4) - t217;
	t214 = cos(qJ(1));
	t208 = pkin(10) + qJ(4);
	t205 = sin(t208);
	t222 = qJD(1) * t205 + qJD(5);
	t237 = t222 * t214;
	t236 = t232 * t205;
	t206 = cos(t208);
	t234 = t232 * t206 - pkin(3) * sin(pkin(10)) - pkin(4) * t205 - qJ(2);
	t212 = sin(qJ(1));
	t227 = qJD(4) * t214;
	t233 = -t206 * t227 + t222 * t212;
	t231 = -pkin(1) - pkin(7) - qJ(3);
	t230 = qJD(1) * t212;
	t229 = qJD(4) * t212;
	t228 = qJD(4) * t213;
	t226 = qJD(5) * t206;
	t225 = qJD(1) * t232;
	t223 = -qJD(5) * t205 - qJD(1);
	t220 = r_i_i_C(1) * t213 - r_i_i_C(2) * t211 + pkin(4);
	t219 = t223 * t214;
	t216 = qJD(1) * t220;
	t215 = qJD(2) + (pkin(4) * t206 + t236) * qJD(4);
	t207 = qJD(1) * t214;
	t204 = t213 * t237 + (t206 * t228 + t223 * t211) * t212;
	t203 = t223 * t213 * t212 + (-t206 * t229 - t237) * t211;
	t202 = t211 * t219 - t233 * t213;
	t201 = t233 * t211 + t213 * t219;
	t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) - t212 * qJD(3) + t215 * t214 + (t234 * t212 + t231 * t214) * qJD(1), t207, -t230, (t214 * t225 - t220 * t229) * t205 + (t238 * t212 + t214 * t216) * t206, t203 * r_i_i_C(1) - t204 * r_i_i_C(2), 0; t204 * r_i_i_C(1) + t203 * r_i_i_C(2) + t214 * qJD(3) + t215 * t212 + (t231 * t212 - t234 * t214) * qJD(1), t230, t207, (t212 * t225 + t220 * t227) * t205 + (t212 * t216 - t238 * t214) * t206, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2), 0; 0, 0, 0, t205 * t217 + (-t220 * t206 - t236) * qJD(4), (t205 * t228 + t211 * t226) * r_i_i_C(2) + (qJD(4) * t205 * t211 - t213 * t226) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:30
	% EndTime: 2019-10-10 00:13:30
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (372->58), mult. (438->90), div. (0->0), fcn. (346->9), ass. (0->52)
	t244 = qJD(5) + qJD(6);
	t248 = sin(qJ(5));
	t286 = pkin(5) * t248;
	t271 = qJD(5) * t286;
	t245 = qJ(5) + qJ(6);
	t241 = sin(t245);
	t242 = cos(t245);
	t294 = r_i_i_C(1) * t241 + r_i_i_C(2) * t242;
	t254 = t294 * t244 + t271;
	t281 = r_i_i_C(3) + pkin(9) + pkin(8);
	t295 = -t281 * qJD(4) + t254;
	t293 = -r_i_i_C(1) * t242 + r_i_i_C(2) * t241;
	t243 = pkin(10) + qJ(4);
	t238 = sin(t243);
	t250 = cos(qJ(5));
	t292 = t250 * (qJD(5) * t238 + qJD(1));
	t291 = t281 * t238;
	t237 = t250 * pkin(5) + pkin(4);
	t239 = cos(t243);
	t289 = t281 * t239 - pkin(3) * sin(pkin(10)) - t237 * t238 - qJ(2);
	t249 = sin(qJ(1));
	t277 = qJD(1) * t238;
	t266 = t244 + t277;
	t251 = cos(qJ(1));
	t273 = qJD(4) * t251;
	t269 = t239 * t273;
	t288 = t266 * t249 - t269;
	t274 = qJD(4) * t249;
	t270 = t239 * t274;
	t287 = t266 * t251 + t270;
	t267 = -t238 * t244 - qJD(1);
	t259 = t267 * t251;
	t230 = t288 * t241 + t242 * t259;
	t231 = t241 * t259 - t288 * t242;
	t279 = -t230 * r_i_i_C(1) + t231 * r_i_i_C(2);
	t260 = t267 * t249;
	t232 = -t287 * t241 + t242 * t260;
	t233 = t241 * t260 + t287 * t242;
	t278 = t232 * r_i_i_C(1) - t233 * r_i_i_C(2);
	t276 = qJD(1) * t249;
	t275 = qJD(4) * t238;
	t272 = qJD(5) * t250;
	t268 = qJD(1) * t281;
	t264 = -qJD(5) - t277;
	t263 = -pkin(1) - pkin(7) - qJ(3) - t286;
	t262 = pkin(5) * t272 + qJD(3);
	t258 = t237 - t293;
	t256 = qJD(1) * t258;
	t255 = t293 * t239 * t244 + t294 * t275;
	t253 = -t238 * t271 + qJD(2) + (t237 * t239 + t291) * qJD(4);
	t240 = qJD(1) * t251;
	t1 = [t231 * r_i_i_C(1) + t230 * r_i_i_C(2) - t262 * t249 + t253 * t251 + (t289 * t249 + t263 * t251) * qJD(1), t240, -t276, (t251 * t268 - t258 * t274) * t238 + (-t295 * t249 + t251 * t256) * t239, (-t249 * t292 + (t264 * t251 - t270) * t248) * pkin(5) + t278, t278; t233 * r_i_i_C(1) + t232 * r_i_i_C(2) + t262 * t251 + t253 * t249 + (t263 * t249 - t289 * t251) * qJD(1), t276, t240, (t249 * t268 + t258 * t273) * t238 + (t249 * t256 + t295 * t251) * t239, (t251 * t292 + (t264 * t249 + t269) * t248) * pkin(5) + t279, t279; 0, 0, 0, t254 * t238 + (-t258 * t239 - t291) * qJD(4), (-t239 * t272 + t248 * t275) * pkin(5) + t255, t255;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end