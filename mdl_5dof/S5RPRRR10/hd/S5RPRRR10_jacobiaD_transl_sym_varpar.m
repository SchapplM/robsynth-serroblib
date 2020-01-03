% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRR10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:11:37
	% EndTime: 2019-12-31 19:11:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:11:37
	% EndTime: 2019-12-31 19:11:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:11:37
	% EndTime: 2019-12-31 19:11:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(9)) + r_i_i_C(2) * sin(pkin(9)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:11:37
	% EndTime: 2019-12-31 19:11:37
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (43->19), mult. (70->32), div. (0->0), fcn. (46->5), ass. (0->14)
	t35 = r_i_i_C(3) + pkin(6) + qJ(2);
	t26 = sin(qJ(1));
	t34 = qJD(1) * t26;
	t27 = cos(qJ(1));
	t33 = qJD(1) * t27;
	t32 = qJD(3) * t26;
	t31 = qJD(3) * t27;
	t24 = pkin(9) + qJ(3);
	t22 = sin(t24);
	t23 = cos(t24);
	t30 = r_i_i_C(1) * t22 + r_i_i_C(2) * t23;
	t29 = -r_i_i_C(1) * t23 + r_i_i_C(2) * t22 - cos(pkin(9)) * pkin(2) - pkin(1);
	t28 = t30 * qJD(3);
	t1 = [t27 * qJD(2) + t30 * t32 + (-t35 * t26 + t29 * t27) * qJD(1), t33, (t22 * t31 + t23 * t34) * r_i_i_C(2) + (t22 * t34 - t23 * t31) * r_i_i_C(1), 0, 0; t26 * qJD(2) - t27 * t28 + (t29 * t26 + t35 * t27) * qJD(1), t34, (t22 * t32 - t23 * t33) * r_i_i_C(2) + (-t22 * t33 - t23 * t32) * r_i_i_C(1), 0, 0; 0, 0, -t28, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:11:38
	% EndTime: 2019-12-31 19:11:38
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (160->43), mult. (276->76), div. (0->0), fcn. (217->7), ass. (0->33)
	t209 = pkin(9) + qJ(3);
	t207 = sin(t209);
	t208 = cos(t209);
	t234 = pkin(7) + r_i_i_C(3);
	t224 = t234 * t208;
	t235 = -pkin(3) * t207 + t224;
	t213 = cos(qJ(4));
	t214 = cos(qJ(1));
	t232 = t213 * t214;
	t212 = sin(qJ(1));
	t231 = qJD(1) * t212;
	t230 = qJD(1) * t214;
	t229 = qJD(3) * t212;
	t228 = qJD(3) * t213;
	t227 = qJD(3) * t214;
	t226 = qJD(4) * t207;
	t225 = qJD(4) * t208;
	t223 = -qJD(1) + t225;
	t222 = qJD(1) * t208 - qJD(4);
	t211 = sin(qJ(4));
	t221 = r_i_i_C(1) * t211 + r_i_i_C(2) * t213;
	t220 = r_i_i_C(1) * t213 - r_i_i_C(2) * t211 + pkin(3);
	t219 = t223 * t211;
	t218 = -pkin(3) * t208 - t234 * t207 - cos(pkin(9)) * pkin(2) - pkin(1);
	t217 = qJD(3) * t220;
	t216 = t207 * t227 + t222 * t212;
	t215 = -t234 * qJD(3) + t221 * qJD(4);
	t210 = -pkin(6) - qJ(2);
	t205 = -t222 * t232 + (t207 * t228 + t219) * t212;
	t204 = t223 * t213 * t212 + (-t207 * t229 + t222 * t214) * t211;
	t203 = t216 * t213 + t214 * t219;
	t202 = t216 * t211 - t223 * t232;
	t1 = [t205 * r_i_i_C(1) + t204 * r_i_i_C(2) + t214 * qJD(2) - t235 * t229 + (t210 * t212 + t218 * t214) * qJD(1), t230, (-t214 * t217 - t234 * t231) * t208 + (t215 * t214 + t220 * t231) * t207, t202 * r_i_i_C(1) + t203 * r_i_i_C(2), 0; -t203 * r_i_i_C(1) + t202 * r_i_i_C(2) + t212 * qJD(2) + t235 * t227 + (-t210 * t214 + t218 * t212) * qJD(1), t231, (-t212 * t217 + t234 * t230) * t208 + (t215 * t212 - t220 * t230) * t207, -t204 * r_i_i_C(1) + t205 * r_i_i_C(2), 0; 0, 0, -t221 * t225 + (-t220 * t207 + t224) * qJD(3), (-t208 * t228 + t211 * t226) * r_i_i_C(2) + (-qJD(3) * t208 * t211 - t213 * t226) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:11:38
	% EndTime: 2019-12-31 19:11:38
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (367->60), mult. (426->91), div. (0->0), fcn. (338->9), ass. (0->53)
	t251 = cos(qJ(4));
	t240 = t251 * pkin(4) + pkin(3);
	t245 = pkin(9) + qJ(3);
	t241 = sin(t245);
	t242 = cos(t245);
	t283 = r_i_i_C(3) + pkin(8) + pkin(7);
	t268 = t283 * t242;
	t249 = sin(qJ(4));
	t282 = pkin(4) * qJD(4);
	t272 = t249 * t282;
	t291 = (-t240 * t241 + t268) * qJD(3) - t242 * t272;
	t252 = cos(qJ(1));
	t246 = qJD(4) + qJD(5);
	t266 = t242 * t246 - qJD(1);
	t289 = t252 * t266;
	t247 = qJ(4) + qJ(5);
	t243 = sin(t247);
	t244 = cos(t247);
	t284 = r_i_i_C(2) * t244;
	t261 = r_i_i_C(1) * t243 + t284;
	t288 = -t261 * t246 - t272;
	t277 = qJD(1) * t242;
	t265 = -t246 + t277;
	t250 = sin(qJ(1));
	t274 = qJD(3) * t241;
	t270 = t250 * t274;
	t287 = t265 * t252 - t270;
	t286 = pkin(4) * t249;
	t285 = r_i_i_C(2) * t243;
	t280 = t244 * t246;
	t269 = t252 * t274;
	t255 = t265 * t250 + t269;
	t234 = t255 * t243 - t244 * t289;
	t235 = t243 * t289 + t255 * t244;
	t279 = t234 * r_i_i_C(1) + t235 * r_i_i_C(2);
	t260 = t266 * t250;
	t236 = t287 * t243 + t244 * t260;
	t237 = t243 * t260 - t287 * t244;
	t278 = -t236 * r_i_i_C(1) + t237 * r_i_i_C(2);
	t276 = qJD(1) * t250;
	t275 = qJD(1) * t252;
	t273 = qJD(3) * t242;
	t271 = t251 * t282;
	t267 = pkin(6) + qJ(2) + t286;
	t263 = -qJD(4) + t277;
	t262 = qJD(2) + t271;
	t259 = (-qJD(4) * t242 + qJD(1)) * t251;
	t258 = r_i_i_C(1) * t244 + t240 - t285;
	t257 = qJD(3) * t258;
	t256 = -t240 * t242 - t283 * t241 - cos(pkin(9)) * pkin(2) - pkin(1);
	t254 = -t283 * qJD(3) - t288;
	t238 = t241 * t246 * t285;
	t1 = [t237 * r_i_i_C(1) + t236 * r_i_i_C(2) + t262 * t252 - t291 * t250 + (-t267 * t250 + t256 * t252) * qJD(1), t275, (-t252 * t257 - t283 * t276) * t242 + (t254 * t252 + t258 * t276) * t241, (t252 * t259 + (t263 * t250 + t269) * t249) * pkin(4) + t279, t279; -t235 * r_i_i_C(1) + t234 * r_i_i_C(2) + t262 * t250 + t291 * t252 + (t256 * t250 + t267 * t252) * qJD(1), t276, (-t250 * t257 + t283 * t275) * t242 + (t254 * t250 - t258 * t275) * t241, (t250 * t259 + (-t263 * t252 + t270) * t249) * pkin(4) + t278, t278; 0, 0, t288 * t242 + (-t258 * t241 + t268) * qJD(3), t238 + (-r_i_i_C(1) * t280 - t271) * t241 + (-t261 - t286) * t273, -t273 * t284 + t238 + (-t241 * t280 - t243 * t273) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end