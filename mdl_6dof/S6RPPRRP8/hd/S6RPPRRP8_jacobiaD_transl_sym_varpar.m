% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:59
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP8_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP8_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP8_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:41
	% EndTime: 2019-10-09 23:59:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:41
	% EndTime: 2019-10-09 23:59:41
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
	% StartTime: 2019-10-09 23:59:41
	% EndTime: 2019-10-09 23:59:41
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
	% StartTime: 2019-10-09 23:59:41
	% EndTime: 2019-10-09 23:59:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (15->10), mult. (36->14), div. (0->0), fcn. (24->4), ass. (0->7)
	t14 = sin(qJ(1));
	t18 = qJD(1) * t14;
	t17 = -pkin(1) - r_i_i_C(3) - qJ(3);
	t16 = r_i_i_C(1) * sin(pkin(9)) + r_i_i_C(2) * cos(pkin(9)) + qJ(2);
	t15 = cos(qJ(1));
	t11 = qJD(1) * t15;
	t1 = [t15 * qJD(2) - t14 * qJD(3) + (-t16 * t14 + t17 * t15) * qJD(1), t11, -t18, 0, 0, 0; t14 * qJD(2) + t15 * qJD(3) + (t17 * t14 + t16 * t15) * qJD(1), t18, t11, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:41
	% EndTime: 2019-10-09 23:59:41
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (48->21), mult. (82->32), div. (0->0), fcn. (54->5), ass. (0->14)
	t24 = pkin(9) + qJ(4);
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
	t30 = pkin(3) * sin(pkin(9)) + r_i_i_C(1) * t21 + r_i_i_C(2) * t22 + qJ(2);
	t29 = qJD(2) + t36;
	t1 = [-t27 * qJD(3) + t29 * t28 + (-t30 * t27 + t32 * t28) * qJD(1), t23, -t35, (-t21 * t23 - t22 * t34) * r_i_i_C(2) + (-t21 * t34 + t22 * t23) * r_i_i_C(1), 0, 0; t28 * qJD(3) + t29 * t27 + (t32 * t27 + t30 * t28) * qJD(1), t35, t23, (-t21 * t35 + t22 * t33) * r_i_i_C(2) + (t21 * t33 + t22 * t35) * r_i_i_C(1), 0, 0; 0, 0, 0, -t36, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:59:42
	% EndTime: 2019-10-09 23:59:42
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (165->45), mult. (288->78), div. (0->0), fcn. (225->7), ass. (0->33)
	t211 = sin(qJ(5));
	t213 = cos(qJ(5));
	t217 = (r_i_i_C(1) * t211 + r_i_i_C(2) * t213) * qJD(5);
	t232 = pkin(8) + r_i_i_C(3);
	t238 = t232 * qJD(4) - t217;
	t214 = cos(qJ(1));
	t208 = pkin(9) + qJ(4);
	t205 = sin(t208);
	t222 = qJD(1) * t205 + qJD(5);
	t237 = t222 * t214;
	t236 = t232 * t205;
	t206 = cos(t208);
	t234 = t232 * t206 - pkin(3) * sin(pkin(9)) - pkin(4) * t205 - qJ(2);
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
	% StartTime: 2019-10-09 23:59:42
	% EndTime: 2019-10-09 23:59:43
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (306->57), mult. (528->87), div. (0->0), fcn. (441->7), ass. (0->44)
	t248 = sin(qJ(5));
	t250 = cos(qJ(5));
	t278 = r_i_i_C(3) + qJ(6);
	t281 = pkin(5) + r_i_i_C(1);
	t258 = -t278 * t248 - t281 * t250;
	t256 = -pkin(4) + t258;
	t245 = pkin(9) + qJ(4);
	t242 = sin(t245);
	t280 = pkin(8) + r_i_i_C(2);
	t284 = t280 * t242;
	t243 = cos(t245);
	t283 = t280 * t243 - pkin(3) * sin(pkin(9)) - pkin(4) * t242 - qJ(2);
	t257 = t281 * t248 - t278 * t250;
	t269 = qJD(6) * t248;
	t253 = t257 * qJD(5) - t269;
	t282 = -t280 * qJD(4) + t253;
	t279 = -pkin(1) - pkin(7) - qJ(3);
	t277 = t243 * t250;
	t249 = sin(qJ(1));
	t276 = t249 * t248;
	t275 = t249 * t250;
	t251 = cos(qJ(1));
	t274 = t251 * t248;
	t273 = qJD(1) * t249;
	t244 = qJD(1) * t251;
	t272 = qJD(4) * t242;
	t271 = qJD(4) * t249;
	t270 = qJD(4) * t251;
	t268 = t250 * qJD(6);
	t267 = t251 * t242 * t250;
	t266 = qJD(1) * t280;
	t265 = t243 * t270;
	t264 = qJD(5) * t242 + qJD(1);
	t263 = qJD(1) * t242 + qJD(5);
	t262 = -qJD(3) + t268;
	t261 = t263 * t251;
	t259 = t242 * t275 + t274;
	t254 = qJD(1) * t256;
	t252 = t242 * t269 + qJD(2) + (pkin(4) * t243 + t284) * qJD(4);
	t237 = t250 * t261 + (qJD(4) * t277 - t264 * t248) * t249;
	t236 = t264 * t275 + (t243 * t271 + t261) * t248;
	t235 = -t250 * t265 + (t242 * t274 + t275) * qJD(5) + t259 * qJD(1);
	t234 = -qJD(5) * t267 - t250 * t244 - t248 * t265 + t263 * t276;
	t1 = [t262 * t249 - t281 * t235 - t278 * t234 + t252 * t251 + (t283 * t249 + t279 * t251) * qJD(1), t244, -t273, (t251 * t266 + t256 * t271) * t242 + (-t282 * t249 - t251 * t254) * t243, t259 * qJD(6) - t281 * t236 + t278 * t237, t236; -t262 * t251 + t281 * t237 + t278 * t236 + t252 * t249 + (t279 * t249 - t283 * t251) * qJD(1), t273, t244, (t249 * t266 - t256 * t270) * t242 + (-t249 * t254 + t282 * t251) * t243, -(t267 - t276) * qJD(6) + t278 * t235 - t281 * t234, t234; 0, 0, 0, t253 * t242 + (t256 * t243 - t284) * qJD(4), t257 * t272 + (t258 * qJD(5) + t268) * t243, qJD(5) * t277 - t248 * t272;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end