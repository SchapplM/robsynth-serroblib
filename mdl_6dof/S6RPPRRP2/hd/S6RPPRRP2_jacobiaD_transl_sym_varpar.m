% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
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
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(9);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (26->10), mult. (32->14), div. (0->0), fcn. (20->6), ass. (0->6)
	t19 = r_i_i_C(3) + qJ(3);
	t18 = -r_i_i_C(1) * cos(pkin(10)) + r_i_i_C(2) * sin(pkin(10)) - pkin(2);
	t15 = qJ(1) + pkin(9);
	t14 = cos(t15);
	t13 = sin(t15);
	t1 = [t14 * qJD(3) + (-cos(qJ(1)) * pkin(1) - t19 * t13 + t18 * t14) * qJD(1), 0, qJD(1) * t14, 0, 0, 0; t13 * qJD(3) + (-sin(qJ(1)) * pkin(1) + t19 * t14 + t18 * t13) * qJD(1), 0, qJD(1) * t13, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (71->22), mult. (74->34), div. (0->0), fcn. (48->7), ass. (0->15)
	t40 = r_i_i_C(3) + pkin(7) + qJ(3);
	t31 = qJ(1) + pkin(9);
	t27 = sin(t31);
	t39 = qJD(1) * t27;
	t29 = cos(t31);
	t38 = qJD(1) * t29;
	t37 = qJD(4) * t27;
	t36 = qJD(4) * t29;
	t30 = pkin(10) + qJ(4);
	t26 = sin(t30);
	t28 = cos(t30);
	t35 = r_i_i_C(1) * t26 + r_i_i_C(2) * t28;
	t34 = -r_i_i_C(1) * t28 + r_i_i_C(2) * t26 - cos(pkin(10)) * pkin(3) - pkin(2);
	t33 = t35 * qJD(4);
	t1 = [t29 * qJD(3) + t35 * t37 + (-cos(qJ(1)) * pkin(1) - t40 * t27 + t34 * t29) * qJD(1), 0, t38, (t26 * t36 + t28 * t39) * r_i_i_C(2) + (t26 * t39 - t28 * t36) * r_i_i_C(1), 0, 0; t27 * qJD(3) - t29 * t33 + (-sin(qJ(1)) * pkin(1) + t40 * t29 + t34 * t27) * qJD(1), 0, t39, (t26 * t37 - t28 * t38) * r_i_i_C(2) + (-t26 * t38 - t28 * t37) * r_i_i_C(1), 0, 0; 0, 0, 0, -t33, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:31
	% EndTime: 2019-10-09 23:49:31
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (246->45), mult. (280->75), div. (0->0), fcn. (219->9), ass. (0->34)
	t217 = pkin(10) + qJ(4);
	t213 = sin(t217);
	t215 = cos(t217);
	t240 = pkin(8) + r_i_i_C(3);
	t232 = t240 * t215;
	t244 = (-pkin(4) * t213 + t232) * qJD(4);
	t221 = cos(qJ(5));
	t229 = qJD(1) * t215 - qJD(5);
	t242 = t221 * t229;
	t233 = qJD(5) * t215;
	t230 = -qJD(1) + t233;
	t220 = sin(qJ(5));
	t236 = qJD(4) * t220;
	t241 = -t213 * t236 + t230 * t221;
	t218 = qJ(1) + pkin(9);
	t214 = sin(t218);
	t238 = qJD(1) * t214;
	t216 = cos(t218);
	t237 = qJD(1) * t216;
	t235 = qJD(4) * t221;
	t234 = qJD(5) * t213;
	t228 = r_i_i_C(1) * t220 + r_i_i_C(2) * t221;
	t227 = r_i_i_C(1) * t221 - r_i_i_C(2) * t220 + pkin(4);
	t226 = t229 * t220;
	t225 = -pkin(4) * t215 - t240 * t213 - cos(pkin(10)) * pkin(3) - pkin(2);
	t224 = qJD(4) * t227;
	t223 = t213 * t235 + t230 * t220;
	t222 = -t240 * qJD(4) + t228 * qJD(5);
	t219 = -pkin(7) - qJ(3);
	t211 = t223 * t214 - t216 * t242;
	t210 = t241 * t214 + t216 * t226;
	t209 = t214 * t242 + t223 * t216;
	t208 = t214 * t226 - t241 * t216;
	t1 = [t211 * r_i_i_C(1) + t210 * r_i_i_C(2) + t216 * qJD(3) - t214 * t244 + (-cos(qJ(1)) * pkin(1) + t214 * t219 + t225 * t216) * qJD(1), 0, t237, (-t216 * t224 - t240 * t238) * t215 + (t222 * t216 + t227 * t238) * t213, t208 * r_i_i_C(1) + t209 * r_i_i_C(2), 0; -t209 * r_i_i_C(1) + t208 * r_i_i_C(2) + t214 * qJD(3) + t216 * t244 + (-sin(qJ(1)) * pkin(1) - t216 * t219 + t225 * t214) * qJD(1), 0, t238, (-t214 * t224 + t240 * t237) * t215 + (t222 * t214 - t227 * t237) * t213, -t210 * r_i_i_C(1) + t211 * r_i_i_C(2), 0; 0, 0, 0, -t228 * t233 + (-t227 * t213 + t232) * qJD(4), (-t215 * t235 + t220 * t234) * r_i_i_C(2) + (-t215 * t236 - t221 * t234) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:32
	% EndTime: 2019-10-09 23:49:32
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (459->61), mult. (520->91), div. (0->0), fcn. (435->9), ass. (0->45)
	t269 = pkin(10) + qJ(4);
	t265 = sin(t269);
	t267 = cos(t269);
	t304 = pkin(8) + r_i_i_C(2);
	t288 = t304 * t267;
	t306 = -pkin(4) * t265 + t288;
	t272 = sin(qJ(5));
	t273 = cos(qJ(5));
	t301 = r_i_i_C(3) + qJ(6);
	t303 = -r_i_i_C(1) - pkin(5);
	t280 = -t301 * t272 + t303 * t273;
	t276 = -pkin(4) + t280;
	t279 = t303 * t272 + t301 * t273;
	t305 = t279 * qJD(5) + qJD(6) * t272;
	t270 = qJ(1) + pkin(9);
	t266 = sin(t270);
	t300 = t266 * t272;
	t299 = t266 * t273;
	t268 = cos(t270);
	t298 = t268 * t272;
	t297 = t268 * t273;
	t296 = qJD(1) * t266;
	t295 = qJD(1) * t268;
	t294 = qJD(4) * t266;
	t293 = qJD(4) * t267;
	t292 = qJD(4) * t268;
	t291 = qJD(5) * t272;
	t290 = qJD(5) * t273;
	t287 = t265 * t294;
	t286 = t266 * t291;
	t285 = t265 * t292;
	t284 = t268 * t290;
	t283 = t267 * t297 + t300;
	t282 = t267 * t300 + t297;
	t281 = -pkin(4) * t267 - t304 * t265 - cos(pkin(10)) * pkin(3) - pkin(2);
	t278 = t266 * t290 + t272 * t295;
	t277 = t268 * t291 + t273 * t296;
	t275 = qJD(4) * t276;
	t274 = -t304 * qJD(4) - t305;
	t271 = -pkin(7) - qJ(3);
	t253 = t283 * qJD(1) - t267 * t286 - t273 * t287 - t284;
	t252 = t278 * t267 - t272 * t287 - t277;
	t251 = t277 * t267 + t273 * t285 - t278;
	t250 = t282 * qJD(1) - t267 * t284 + t272 * t285 - t286;
	t1 = [-t282 * qJD(6) + t268 * qJD(3) + t303 * t253 - t301 * t252 - t306 * t294 + (-cos(qJ(1)) * pkin(1) + t266 * t271 + t281 * t268) * qJD(1), 0, t295, (t268 * t275 - t304 * t296) * t267 + (t274 * t268 - t276 * t296) * t265, t283 * qJD(6) - t303 * t250 - t301 * t251, -t250; -(-t267 * t298 + t299) * qJD(6) + t266 * qJD(3) + t303 * t251 - t301 * t250 + t306 * t292 + (-sin(qJ(1)) * pkin(1) - t268 * t271 + t281 * t266) * qJD(1), 0, t296, (t266 * t275 + t304 * t295) * t267 + (t274 * t266 + t276 * t295) * t265, -(-t267 * t299 + t298) * qJD(6) + t301 * t253 + t303 * t252, t252; 0, 0, 0, t305 * t267 + (t276 * t265 + t288) * qJD(4), t279 * t293 + (t280 * qJD(5) + t273 * qJD(6)) * t265, t265 * t290 + t272 * t293;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end