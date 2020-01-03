% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPP7
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPP7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPP7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_jacobiaD_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:06:38
	% EndTime: 2019-12-31 21:06:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:06:38
	% EndTime: 2019-12-31 21:06:39
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
	% StartTime: 2019-12-31 21:06:38
	% EndTime: 2019-12-31 21:06:38
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(6) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t18 * t28 + t20 * t22) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0; -t20 * t21 + (t18 * t22 + t20 * t28) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0; 0, -t21, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:06:39
	% EndTime: 2019-12-31 21:06:40
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (83->38), mult. (270->73), div. (0->0), fcn. (211->6), ass. (0->31)
	t199 = sin(qJ(2));
	t202 = cos(qJ(2));
	t223 = pkin(7) + r_i_i_C(3);
	t213 = t223 * t202;
	t224 = -pkin(2) * t199 + t213;
	t201 = cos(qJ(3));
	t203 = cos(qJ(1));
	t221 = t201 * t203;
	t200 = sin(qJ(1));
	t220 = qJD(1) * t200;
	t219 = qJD(1) * t203;
	t218 = qJD(2) * t200;
	t217 = qJD(2) * t202;
	t216 = qJD(2) * t203;
	t215 = qJD(3) * t199;
	t214 = qJD(3) * t202;
	t212 = -qJD(1) + t214;
	t211 = qJD(1) * t202 - qJD(3);
	t198 = sin(qJ(3));
	t210 = r_i_i_C(1) * t198 + r_i_i_C(2) * t201;
	t209 = r_i_i_C(1) * t201 - r_i_i_C(2) * t198 + pkin(2);
	t208 = t212 * t198;
	t207 = -pkin(2) * t202 - t223 * t199 - pkin(1);
	t206 = qJD(2) * t209;
	t205 = t199 * t216 + t211 * t200;
	t204 = -t223 * qJD(2) + t210 * qJD(3);
	t197 = -t211 * t221 + (qJD(2) * t199 * t201 + t208) * t200;
	t196 = t212 * t201 * t200 + (-t199 * t218 + t211 * t203) * t198;
	t195 = t205 * t201 + t203 * t208;
	t194 = t205 * t198 - t212 * t221;
	t1 = [t197 * r_i_i_C(1) + t196 * r_i_i_C(2) - t224 * t218 + (-pkin(6) * t200 + t207 * t203) * qJD(1), (-t203 * t206 - t223 * t220) * t202 + (t204 * t203 + t209 * t220) * t199, t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0, 0; -t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t224 * t216 + (pkin(6) * t203 + t207 * t200) * qJD(1), (-t200 * t206 + t223 * t219) * t202 + (t204 * t200 - t209 * t219) * t199, -t196 * r_i_i_C(1) + t197 * r_i_i_C(2), 0, 0; 0, -t210 * t214 + (-t209 * t199 + t213) * qJD(2), (t198 * t215 - t201 * t217) * r_i_i_C(2) + (-t198 * t217 - t201 * t215) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:06:40
	% EndTime: 2019-12-31 21:06:40
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (164->54), mult. (510->89), div. (0->0), fcn. (427->6), ass. (0->41)
	t248 = sin(qJ(3));
	t251 = cos(qJ(3));
	t280 = r_i_i_C(3) + qJ(4);
	t283 = pkin(3) + r_i_i_C(1);
	t284 = t283 * t248 - t280 * t251;
	t286 = -t284 * qJD(3) + qJD(4) * t248;
	t249 = sin(qJ(2));
	t252 = cos(qJ(2));
	t282 = pkin(7) + r_i_i_C(2);
	t268 = t282 * t252;
	t285 = -pkin(2) * t249 + t268;
	t259 = -t280 * t248 - t283 * t251;
	t256 = -pkin(2) + t259;
	t250 = sin(qJ(1));
	t279 = t250 * t248;
	t278 = t250 * t252;
	t253 = cos(qJ(1));
	t277 = t252 * t253;
	t276 = qJD(1) * t250;
	t275 = qJD(1) * t253;
	t274 = qJD(2) * t250;
	t273 = qJD(2) * t252;
	t272 = qJD(2) * t253;
	t271 = qJD(3) * t251;
	t270 = qJD(3) * t253;
	t267 = t249 * t274;
	t266 = t249 * t272;
	t265 = qJD(3) * t279;
	t264 = t248 * t270;
	t263 = t251 * t270;
	t262 = t251 * t277 + t279;
	t261 = t248 * t278 + t251 * t253;
	t260 = -pkin(2) * t252 - t282 * t249 - pkin(1);
	t257 = t248 * t275 + t250 * t271;
	t255 = qJD(2) * t256;
	t254 = -t282 * qJD(2) - t286;
	t237 = t262 * qJD(1) - t251 * t267 - t252 * t265 - t263;
	t236 = -t248 * t267 - t251 * t276 + t257 * t252 - t264;
	t235 = t252 * t264 + (t252 * t276 + t266) * t251 - t257;
	t234 = t261 * qJD(1) + t248 * t266 - t252 * t263 - t265;
	t1 = [-t261 * qJD(4) - t283 * t237 - t280 * t236 - t285 * t274 + (-t250 * pkin(6) + t260 * t253) * qJD(1), (t253 * t255 - t282 * t276) * t252 + (t254 * t253 - t256 * t276) * t249, t262 * qJD(4) + t283 * t234 - t280 * t235, -t234, 0; -(-t248 * t277 + t250 * t251) * qJD(4) - t283 * t235 - t280 * t234 + t285 * t272 + (t253 * pkin(6) + t260 * t250) * qJD(1), (t250 * t255 + t282 * t275) * t252 + (t254 * t250 + t256 * t275) * t249, -(t248 * t253 - t251 * t278) * qJD(4) + t280 * t237 - t283 * t236, t236, 0; 0, t286 * t252 + (t256 * t249 + t268) * qJD(2), -t284 * t273 + (t259 * qJD(3) + t251 * qJD(4)) * t249, t248 * t273 + t249 * t271, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:06:39
	% EndTime: 2019-12-31 21:06:40
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (211->57), mult. (646->87), div. (0->0), fcn. (540->6), ass. (0->44)
	t214 = sin(qJ(2));
	t217 = cos(qJ(2));
	t239 = t214 * qJD(5);
	t213 = sin(qJ(3));
	t240 = qJD(4) * t213;
	t236 = pkin(7) - r_i_i_C(3) - qJ(5);
	t254 = t236 * t217;
	t257 = (-pkin(2) * t214 + t254) * qJD(2) + t217 * t240 - t239;
	t216 = cos(qJ(3));
	t237 = r_i_i_C(1) + pkin(4) + pkin(3);
	t251 = r_i_i_C(2) + qJ(4);
	t253 = t237 * t213 - t251 * t216;
	t256 = t253 * qJD(3) - t240;
	t223 = -t251 * t213 - t237 * t216;
	t221 = -pkin(2) + t223;
	t220 = t221 * t214 + t254;
	t215 = sin(qJ(1));
	t250 = t215 * t213;
	t249 = t215 * t217;
	t218 = cos(qJ(1));
	t248 = t216 * t218;
	t247 = qJD(1) * t215;
	t246 = qJD(1) * t218;
	t245 = qJD(2) * t214;
	t244 = qJD(2) * t217;
	t243 = qJD(2) * t218;
	t242 = qJD(3) * t216;
	t241 = qJD(3) * t218;
	t238 = t216 * qJD(4);
	t235 = t215 * t245;
	t234 = t214 * t243;
	t233 = qJD(3) * t250;
	t232 = t213 * t241;
	t231 = t216 * t241;
	t230 = t236 * t214;
	t227 = t217 * t248 + t250;
	t225 = t213 * t246 + t215 * t242;
	t224 = -pkin(2) * t217 - pkin(1) - t230;
	t219 = -qJD(5) * t217 + t256 * t214 + (t221 * t217 - t230) * qJD(2);
	t202 = t227 * qJD(1) - t216 * t235 - t217 * t233 - t231;
	t201 = -t213 * t235 - t216 * t247 + t225 * t217 - t232;
	t200 = t217 * t232 + (t217 * t247 + t234) * t216 - t225;
	t199 = t213 * t234 - t217 * t231 - t233 + (t213 * t249 + t248) * qJD(1);
	t1 = [-t218 * t238 - t251 * t201 - t237 * t202 - t257 * t215 + (-t215 * pkin(6) + t224 * t218) * qJD(1), t219 * t218 - t220 * t247, t227 * qJD(4) + t237 * t199 - t251 * t200, -t199, t214 * t247 - t217 * t243; -t215 * t238 - t251 * t199 - t237 * t200 + t257 * t218 + (t218 * pkin(6) + t224 * t215) * qJD(1), t219 * t215 + t220 * t246, -(t218 * t213 - t216 * t249) * qJD(4) + t251 * t202 - t237 * t201, t201, -t214 * t246 - t215 * t244; 0, t220 * qJD(2) - t256 * t217 - t239, -t253 * t244 + (t223 * qJD(3) + t238) * t214, t213 * t244 + t214 * t242, -t245;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end